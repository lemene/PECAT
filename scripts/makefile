.PHONY: scripts

TARGET_DIR=../build/bin

TARGET_FILES =   ${TARGET_DIR}/FsaUtils.pm 
TARGET_SCRIPTS = ${TARGET_DIR}/pecat.pl ${TARGET_DIR}/fxtools.py \
				 ${TARGET_DIR}/phase_using_whatshap.py ${TARGET_DIR}/misc.py ${TARGET_DIR}/show_hapkmer.R \
				 ${TARGET_DIR}/prjfile.py ${TARGET_DIR}/fxbug.py ${TARGET_DIR}/parallel_medaka.py

scripts: ${TARGET_SCRIPTS} ${TARGET_FILES}


${TARGET_SCRIPTS}:${TARGET_DIR}/% : ./% 
	cp -pf  $< $@
	chmod +x $@

${TARGET_FILES}:${TARGET_DIR}/% : ./% 
	cp -pf  $< $@

clean:
	rm -f ${TARGET_SCRIPTS} ${TARGET_FILES}


