Input Specification
%%%%%%%%%%%%%%%%%%%

ATS input files are xml files, which can be quite annoying for a human
to create from scratch.  xml is great from the code's perspective, but
it can be quite unintuitive to a human.  That said, is it really that
much worse than structured text with meaningful whitespace and other
nightmares from other codebases?  There really isn't a perfect
solution yet.

We have several approaches in the works to simplify working with ATS
input files.  This includes python-based front-ends for creating them,
higher level specs using YAML or other formats, and more.
Contributions in this front would be fantastic!  In the meantime, the
input spec is what it is, making this documentation critical.

The input spec changes slightly from version to version, but an input
file should never stop working due to commits within the same ATS
release version.  Furthermore, whenever possible, we try to include
python scripts for taking an input file that worked in one version and
updating it to the next newer version.  These are available at
`<https://github.com/amanzi/ats/tree/master/tools/input_converters>`_

.. toctree::
   :maxdepth: 1
   :caption: Versions
             
   input_spec/ATSNativeSpec_dev
   input_spec/ATSNativeSpec_1_1
   input_spec/ATSNativeSpec_1_0
   input_spec/ATSNativeSpec_0_86
