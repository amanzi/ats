#
#  Makefile: generating rst from rst.in

SRC = $(shell find . -name \*.in | cut -c 3-)
RST = $(SRC:%.rst.in=%.rst)

# generic rules for targets

all : rst

%.rst:
	./generate_ATSNativeSpec.py $< $@
	grep -q '^$@' .gitignore || echo "$@" >> .gitignore

rst : $(RST)

clean :
	-/bin/rm -f $(RST)
	-/bin/rm -f *~

