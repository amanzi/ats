# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

.PHONY: all sphinx_help worktree html pdf input_spec deploy clean

all: input_spec html


sphinx_help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

worktree:
	cd build && git worktree add html gh-pages

html:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

pdf:
	@$(SPHINXBUILD) -M pdf "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

input_spec:
	$(MAKE) -C source/input_spec clean rst

deploy:
	cd build/html && \
	git add --all && \
	git commit -m "Deployment to gh-pages" && \
	git push origin gh-pages

clean:
	rm -rf build/html/*

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
#%: Makefile
#	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

