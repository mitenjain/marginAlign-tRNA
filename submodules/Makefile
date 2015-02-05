# order is important, libraries first
modules = sonLib jobTree bwa last cactus
.PHONY: all %.all clean %.clean

all : tokyocabinetRule ${modules:%=all.%} 

tokyocabinetRule :
	cd tokyocabinet &&./configure --prefix=$(PWD)/tokyocabinet ${ktlinkingflags} --disable-bzip --with-zlib=$(PWD)/zlib && make && make install
	export PATH=$(PWD)/tokyocabinet/bin:$PATH
	export tokyoCabinetIncl="-I$(PWD)/tokyocabinet/include -DHAVE_TOKYO_CABINET=1"
	export tokyoCabinetLib="-L$(PWD)/tokyocabinet/lib -Wl,-rpath,$(PWD)/tokyocabinet/lib -ltokyocabinet -lz -lbz2 -lpthread -lm"

all.%:
	export PATH=$(PWD)/tokyocabinet/bin:$PATH
	export tokyoCabinetIncl="-I$(PWD)/tokyocabinet/include -DHAVE_TOKYO_CABINET=1"
	export tokyoCabinetLib="-L$(PWD)/tokyocabinet/lib -Wl,-rpath,$(PWD)/tokyocabinet/lib -ltokyocabinet -lz -lbz2 -lpthread -lm"
	cd $* && make

clean:  ${modules:%=clean.%}

clean.%:
	cd $* && make clean
