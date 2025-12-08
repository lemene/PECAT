
BUILD_BIN_DIR = ./build/bin
BUILD_OBJ_DIR = ./build/obj
BUILD_DIR += $(BUILD_BIN_DIR) $(BUILD_OBJ_DIR)



.PHONY: all src thirdparty  scripts clean
all: $(BUILD_DIR) src thirdparty scripts $(bin) $(libfsa)

$(BUILD_DIR) :
	mkdir -p $(BUILD_DIR) 

src: thirdparty
	cd ./src && make -j 4 && cd -

thirdparty:
	cd ./thirdparty && make -j && cd -

scripts:
	cd ./scripts && make && cd -

clean:
	cd ./src && make clean && cd -
	cd ./thirdparty && make clean && cd -
	

