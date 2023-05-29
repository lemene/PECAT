#!/bin/bash

# compile the project
cd src && make && cd ..

date_str=`date +%Y%m%d`
bin_path=build/bin

commit_ts=`git log -1 --format="%ct"`
commit_time=`date -d@$commit_ts +"%Y-%m-%d %H:%M:%S"`
current_time=`date +"%Y-%m-%d %H:%M:%S"`
git_version=`git log -1 --format="%H"`

echo "Version: $git_version" > $bin_path/VERSION
echo "Commit Time: $commit_time" >> $bin_path/VERSION
echo "Build Time: $current_time" >> $bin_path/VERSION

tar --transform 's,^,PECAT/,S' -czf build/pecat_0.0.2_${git_version:0:7}.tar.gz $bin_path 
