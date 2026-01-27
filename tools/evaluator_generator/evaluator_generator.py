import sys,os
from sympy.printing import ccode

_template_directory = os.path.dirname(os.path.abspath(__file__))
_templates = {}

def loadTemplate(tname):
    with open(os.path.join(_template_directory, "templates", tname),'r') as tfid:
        _templates[tname] = tfid.read()[:-1] # python seems to force end of file with newline character even though it is not there?

def render(tname, d):
    try:
        template = _templates[tname]
    except KeyError:
        loadTemplate(tname)
        template = _templates[tname]        

    assert type(d) is dict        
    return template.format(**d)
    
class EvalGen(object):
    def __init__(self, name, namespace, descriptor, my_key=None, expression=None,
                 doc=None, **kwargs):
        self.d = {}
        self.setName(name, **kwargs)
        self.setNamespace(namespace, **kwargs)
        self.d['evalNameString'] = descriptor
        if my_key is None:
            my_key = name
        self.d['myKeyFirst'] = my_key.split("_")[0]
        self.d['myKey'] = my_key
        self.d['myKeyMethod'] = ''.join([word[0].upper()+word[1:] for word in my_key.split("_")])
        self.args = []
        self.vars = []
        self.pars = []
        self.par_names = []
        self.par_defaults = []
        self.expression = expression
        if doc is not None:
            self.d['docDict'] = doc
        else:
            self.d['docDict'] = ""

    def setName(self, name, **kwargs):
        self._name = name
        if type(name) is str:
            self.d['evalName'] = name
            self.d['evalClassName'] = ''.join([word[0].upper()+word[1:] for word in name.split("_")])
            self.d['evalNameCaps'] = name.upper()

        else:
            self.d['evalName'] = '_'.join(name)
            self.d['evalClassName'] = ''.join([n[0].upper() + n[1:] for n in name])
            self.d['evalNameCaps'] = ''.join(name).upper()

        if 'evalClassName' in kwargs.keys():
            self.d['evalClassName'] = kwargs['evalClassName']


    def setNamespace(self, namespace, **kwargs):
        self.d['namespace'] = namespace[0].upper()+namespace[1:]
        self.d['namespaceCaps'] = namespace.upper()

    def addArg(self, argname, varname=None):
        if varname is None:
            varname = argname
        self.args.append(argname)
        self.vars.append(varname)

    def addParam(self, dname, dtype, parname, default=None):
        self.pars.append((dtype, dname+"_"))
        self.par_names.append(parname)
        self.par_defaults.append(default)

    def renderCopyConstructor(self):
        return '\n'.join([render('evaluator_keyCopyConstructor.cc', dict(arg=arg,var=var)) for arg,var in zip(self.args,self.vars)])

    def renderKeyDeclaration(self):
        return '\n'.join([render('evaluator_keyDeclaration.hh', dict(arg=arg,var=var)) for arg,var in zip(self.args, self.vars)])

    def renderKeyInitialize(self):
        dicts = []
        for arg,var in zip(self.args, self.vars):
            dicts.append(dict(arg=arg, argString=arg.replace("_", " "), var=var))
        return '\n\n'.join([render('evaluator_keyInitialize.cc', argdict) for argdict in dicts])

    def renderKeyCompositeVector(self):
        return '\n'.join([render('evaluator_keyCompositeVector.cc', dict(arg=arg,var=var)) for arg,var in zip(self.args,self.vars)])        

    def renderKeyEpetraVector(self):
        return '\n'.join([render('evaluator_keyEpetraVector.cc', dict(arg=arg,var=var)) for arg,var in zip(self.args,self.vars)])        
    def renderKeyEpetraVectorIndented(self):
        return '\n'.join([render('evaluator_keyEpetraVectorIndented.cc', dict(arg=arg,var=var)) for arg,var in zip(self.args,self.vars)])        

    def renderMyMethodArgs(self):
        return ", ".join(["%s_v[0][i]"%var for var in self.vars])

    def renderMyMethodDeclarationArgs(self):
        return ", ".join(["double %s"%var for var in self.vars])

    def renderEvaluateModel(self):
        d = dict()
        d['keyEpetraVectorList'] = self.renderKeyEpetraVector()
        d['myKeyMethod'] = self.d['myKeyMethod']
        d['myMethodArgs'] = self.renderMyMethodArgs()
        return render('evaluator_evaluateModel.cc', d)

    def renderEvaluateDerivs(self):
        wrt_list = []

        def getDict(arg,var):
            d = dict(arg=arg,var=var)
            d['keyEpetraVectorList'] = self.renderKeyEpetraVectorIndented()
            d['myKeyMethod'] = self.d['myKeyMethod']
            d['wrtMethod'] = ''.join([word[0].upper()+word[1:] for word in arg.split("_")])
            d['myMethodArgs'] = self.renderMyMethodArgs()
            return d
        
        if len(self.args) > 0:
            d = getDict(self.args[0], self.vars[0])
            d['if_elseif'] = render('evaluator_ifWRT.cc', d)
            wrt_list.append(render('evaluator_evaluateDerivs.cc', d))

        if len(self.args) > 1:
            for arg,var in zip(self.args[1:],self.vars[1:]):
                d = getDict(arg,var)
                d['if_elseif'] = render('evaluator_elseifWRT.cc', d)
                wrt_list.append(render('evaluator_evaluateDerivs.cc', d))

        wrt_list.append('\n'.join(["  } else {",
                                   "    ASSERT(0);",
                                   "  }"]))

        return "\n\n".join(wrt_list)
                            
    def renderModelMethodDeclaration(self):
        return render('model_declaration.hh', dict(myMethod=self.d['myKeyMethod'],
                                                   myMethodDeclarationArgs=self.d['myMethodDeclarationArgs']))

    def renderModelDerivDeclarations(self):
        return '\n'.join([render('model_declaration.hh',
                                 dict(myMethod="D%sD%s"%(self.d['myKeyMethod'],''.join([word[0].upper()+word[1:] for word in arg.split("_")])),
                                      myMethodDeclarationArgs=self.d['myMethodDeclarationArgs'])) for arg in self.args])

    def renderModelMethodImplementation(self):
        if self.expression is not None:
            implementation = ccode(self.expression)
        else:
            implementation = "ASSERT(False)"
        return render('model_methodImplementation.cc', dict(evalClassName=self.d['evalClassName'],
                                                            myMethod=self.d['myKeyMethod'],
                                                            myMethodDeclarationArgs=self.d['myMethodDeclarationArgs'],
                                                            myMethodImplementation=implementation))

    def renderModelDerivImplementations(self):
        impls = []

        for arg,var in zip(self.args,self.vars):
            if self.expression is not None:
                print("differentiation of", self.expression, "with respect to", var)
                implementation = ccode(self.expression.diff(var))
            else:
                implementation = "ASSERT(False)"
            impls.append(render('model_methodImplementation.cc',
                                dict(evalClassName=self.d['evalClassName'],
                                     myMethod="D%sD%s"%(self.d['myKeyMethod'],''.join([word[0].upper()+word[1:] for word in arg.split("_")])),
                                     myMethodDeclarationArgs=self.d['myMethodDeclarationArgs'],
                                     myMethodImplementation=implementation)))
        return '\n\n'.join(impls)
    
    def renderModelParamDeclarations(self):
        return '\n'.join(['  %s %s;'%p for p in self.pars])

    def renderModelParamInitializations(self):
        p_inits = []
        for p, pname, pdefault in zip(self.pars, self.par_names, self.par_defaults):
            if pdefault is not None:
                p_inits.append('  %s = plist.get<%s>("%s", %s);'%(p[1],p[0],pname, str(pdefault)))
            else:
                p_inits.append('  %s = plist.get<%s>("%s");'%(p[1],p[0],pname))

        return '\n'.join(p_inits)

    def genArgs(self):
        # dependencies
        self.d['keyDeclarationList'] = self.renderKeyDeclaration()
        self.d['keyCopyConstructorList'] = self.renderCopyConstructor()
        print(self.d['keyCopyConstructorList'])
        self.d['keyInitializeList'] = self.renderKeyInitialize()
        self.d['keyCompositeVectorList'] = self.renderKeyCompositeVector()
        self.d['myMethodArgs'] = self.renderMyMethodArgs()
        self.d['myMethodDeclarationArgs'] = self.renderMyMethodDeclarationArgs()
        self.d['evaluateModel'] = self.renderEvaluateModel()
        self.d['evaluateDerivs'] = self.renderEvaluateDerivs()

        self.d['modelMethodDeclaration'] = self.renderModelMethodDeclaration()
        self.d['modelDerivDeclarationList'] = self.renderModelDerivDeclarations()
        self.d['paramDeclarationList'] = self.renderModelParamDeclarations()

        self.d['modelMethodImplementation'] = self.renderModelMethodImplementation()
        self.d['modelDerivImplementationList'] = self.renderModelDerivImplementations()
        self.d['modelInitializeParamsList'] = self.renderModelParamInitializations()

def generate_evaluator(name, namespace, descriptor, my_key, dependencies, parameters, **kwargs):
    """Generates an evaluator whose class is [name]Evaluator and model is [name]Model.

    This evaluator/model live in [namespace] and have key (for use in input spec) [descriptor].
    It evaluates the seconadary variable based on [dependencies], using [parameters]
    in the model.

    Inputs:
      name: name of the class, must meet C++ standards for class names

      namespace: containing context of the class, must meet C++ standards for namespaces

      descriptor: short (several word) string used in the input spec to refer to this evaluator

      my_key: key identifier of the variable to be evaluated, in base form (note this may be
              modified by a domain, etc, but is the base method used in the model to evaluate
              the model), e.g. "density"

      dependencies: list of strings or tuples, strings are the key identifier of the dependency,
                     or a tuple with (key identifier, variable name) where variable name is a
                     shorter version of the identifier used in argument lists, etc (must meet
                     C++ variable standards)

      parameters: list of tuples (name, type, descriptor) of parameters to be used by the model,
                  where name is the (C++ standard) variable name, type is one of {double, int,
                  bool, string, Array(double), ...} and descriptor is short string used in
                  pulling the value from the input spec.

      directory: directory where output files are created

    Outputs:
      writes files: [name]_evaluator.hh
                    [name]_evaluator.cc
                    [name]_evaluator_reg.hh
                    [name]_model.hh
                    [name]_model.cc
    """
    eg = EvalGen(name, namespace, descriptor, my_key, **kwargs)
    for dep in dependencies:
        if type(dep) is str:
            eg.addArg(dep)
        else:
            eg.addArg(*dep)

    for par in parameters:
        eg.addParam(*par)

    eg.genArgs()

    try:
        directory = kwargs['directory']
    except KeyError:
        directory = "."

    files = ["evaluator.hh", "evaluator.cc", "evaluator_reg.hh", "model.hh", "model.cc"]
    for outfile in files:
        with open(os.path.join(directory, "%s_%s"%(name,outfile)), 'w') as fid:
            fid.write(render(outfile, eg.d))


      
        
if __name__ == "__main__":
    eg = EvalGen("iem", "energy", "internal energy", "internal_energy", evalClassName='IEM')
    eg.addArg("temperature", "temp");
    eg.addParam("cv", "double", "specific heat capacity")
    eg.genArgs()

    print (render("evaluator.hh", eg.d))
    print (render("evaluator.cc", eg.d))
    print (render("evaluator_reg.hh", eg.d))
    print (render("model.hh", eg.d))
    print (render("model.cc", eg.d))    
