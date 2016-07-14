#!/bin/bash

for file in test*.py; do
	python "${file}"
	echo "${file} done!"
done
