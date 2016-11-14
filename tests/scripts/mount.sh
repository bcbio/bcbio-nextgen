#source this file.


_BUCKET="testbcbio"
checkmounted(){ grep -s "$1" "/proc/mounts";}


if checkmounted "$_BUCKET"; then
	echo "S3 bucket is mounted."
else 
    echo "Starting syslogd..."
    service rsyslog start
	echo "Mounting the bucket..."
	cmd="goofys --sse $_BUCKET /mnt/$_BUCKET"
	echo $cmd
	$cmd
    sleep 1
	if checkmounted "$_BUCKET"; then
		echo "Successfully mounted S3 bucket."
	else
		echo "Something went wrong with the mount."
		echo "That's all I know:"
		sudo grep "fuse\|goofys" /var/log/syslog | tail
	fi
fi

export BCBIO_WORKDIR=/mnt/$_BUCKET/testworkdir
echo "BCBIO_WORKDIR is now $BCBIO_WORKDIR"

deactivate () {
	echo "BCBIO_WORKDIR is no longer $BCBIO_WORKDIR"
	unset BCBIO_WORKDIR
}
