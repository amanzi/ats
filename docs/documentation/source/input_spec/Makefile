#
#  Makefile: generating rst from rst.in

SRC = $(shell find . -name \*.in)
RST = $(SRC:%.rst.in=%.rst)

# generic rules for targets

all : rst

%.rst:
	./generate_ATSNativeSpec.py $< $@

rst : $(RST)

clean :
	-/bin/rm -f $(RST)
	-/bin/rm -f *~

