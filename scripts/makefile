.PHONY: scripts

TARGET_DIR=../build/bin

TARGET_FILES =   ${TARGET_DIR}/FsaUtils.pm 
TARGET_SCRIPTS = ${TARGET_DIR}/pecat.pl  ${TARGET_DIR}/raec.pl ${TARGET_DIR}/show_hapkmer.R \
	${TARGET_DIR}/parallel_medaka.py ${TARGET_DIR}/test_pol.py

scripts: ${TARGET_SCRIPTS} ${TARGET_FILES}


${TARGET_SCRIPTS}:${TARGET_DIR}/% : ./% 
	cp -pf  $< $@
	chmod +x $@

${TARGET_FILES}:${TARGET_DIR}/% : ./% 
	cp -pf  $< $@

clean:
	rm -f ${TARGET_SCRIPTS} ${TARGET_FILES}


