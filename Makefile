SHELL:=/bin/bash

.DEFAULT_GOAL := all

ROOT_DIR:=$(shell dirname "$(realpath $(firstword $(MAKEFILE_LIST)))")
OUTPUT_DIRECTORY=${ROOT_DIR}/output
SUBMODULES_PATH?=${ROOT_DIR}

.EXPORT_ALL_VARIABLES:
DOCKER_BUILDKIT?=1
DOCKER_CONFIG?=

PROJECT:=optinlc
TAG := $(shell git rev-parse --short HEAD)

.PHONY: show-hash
show-hash:
	@echo "Git Short Hash: $(GIT_SHORT_HASH)"



.PHONY: all
all: build

.PHONY: init_submodules
init_submodules:
ifeq ($(wildcard ${SUBMODULES_PATH}/mathematics_toolbox/*),)
	git submodule update --init mathematics_toolbox
else
	@echo "Submodules already initialized, skipping submodule init."
endif

.PHONY: build_mathematics_toolbox
build_mathematics_toolbox:
	cd "${SUBMODULES_PATH}/mathematics_toolbox" && make build

.PHONY: build
build: clean init_submodules build_mathematics_toolbox _build

.PHONY: _build
_build:
	docker build --network host \
                 --tag ${PROJECT}:${TAG} \
                 --build-arg PROJECT=${PROJECT} .
	docker cp $$(docker create --rm ${PROJECT}:${TAG}):/tmp/OptiNLC/OptiNLC/build "${ROOT_DIR}/OptiNLC"

.PHONY: test
test: build
	mkdir -p ${OUTPUT_DIRECTORY}
	#docker run -it -v ${OUTPUT_DIRECTORY}:/tmp/output ${PROJECT}:${TAG} /bin/bash -c 'mkdir -p /tmp/output/ && cd /tmp/output/'
	docker run -it -v ${OUTPUT_DIRECTORY}:/tmp/output ${PROJECT}:${TAG} /bin/bash -c 'cd /tmp/output && /tmp/OptiNLC/OptiNLC/build/OptiNLC_TestRunner -d yes'
	#docker run --volume ${OUTPUT_DIRECTORY}:/tmp/output ${PROJECT}:${TAG} /bin/bash -c 'touch /tmp/output/hello.txt'

.PHONY: plot
plot: 
	gnuplot eigen_plot.gnuplot

.PHONY: run
run: build
	docker run -it ${PROJECT}:${TAG} /tmp/OptiNLC/OptiNLC/build/OptiNLC

.PHONY: clean
clean:  ## Clean OptiNLC build artifacts 
	rm -rf "${ROOT_DIR}/${PROJECT}/build"
	rm -rf "${OUTPUT_DIRECTORY}"
	docker rm $$(docker ps -a -q --filter "ancestor=${PROJECT}:${TAG}") --force 2> /dev/null || true
	docker rmi $$(docker images -q ${PROJECT}:${TAG}) --force 2> /dev/null || true
	docker rmi --force $$(docker images --filter "dangling=true" -q --no-trunc) 2> /dev/null || true

