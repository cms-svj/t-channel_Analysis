source ~/.bash_profile
storage_dir=$(readlink -f $PWD)
export TCHANNEL_BASE=${storage_dir}
source coffeaenv/bin/activate
