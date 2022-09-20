.PHONY: thirdparty pigz 



TARGET_DIR=../build/bin

# pldg
PLGD_PM = $(wildcard ./plgd/Plgd/*.pm) $(wildcard ./plgd/Plgd/Grid/*.pm)  $(wildcard ./plgd/Plgd/Job/*.pm)
TARGET_PLGD_PM = $(patsubst ./plgd/Plgd/%, ${TARGET_DIR}/Plgd/%, $(PLGD_PM))

thirdparty: pigz pldg libksw2

pldg: ${TARGET_PLGD_PM}
	echo ${TARGET_PLGD_PM}

${TARGET_PLGD_PM}:${TARGET_DIR}/Plgd/% : ./plgd/Plgd/% 
	@if [ ! -e ${TARGET_DIR}/Plgd ] ; then mkdir -p ${TARGET_DIR}/Plgd ; fi
	@if [ ! -e ${TARGET_DIR}/Plgd/Grid ] ; then mkdir -p ${TARGET_DIR}/Plgd/Grid ; fi
	@if [ ! -e ${TARGET_DIR}/Plgd/Job ] ; then mkdir -p ${TARGET_DIR}/Plgd/Job ; fi
	cp -pf  $< $@

pigz: ${TARGET_DIR}/pigz


${TARGET_DIR}/pigz: ${TARGET_DIR}/../pigz-2.4/pigz
	cp -pf $< $@

${TARGET_DIR}/../pigz-2.4/pigz: ./pigz/pigz-2.4.tar.gz
	tar -xzf $< -C ${TARGET_DIR}/.. || exit 255
	cd ${TARGET_DIR}/../pigz-2.4 && make



libksw2: ${TARGET_DIR}/libksw2.a


${TARGET_DIR}/libksw2.a: ./ksw2/libksw2.a
	cp -pf $< $@

./ksw2/libksw2.a:
	cd ksw2 && make
	ar -r $@ ./ksw2/*.o
	
clean:
	rm -rf ${TARGET_DIR}/Plgd 
	rm -f ${TARGET_DIR}/pigz
	rm -rf ${TARGET_DIR}/../pigz-2.4
