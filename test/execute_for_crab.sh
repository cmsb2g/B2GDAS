#!/usr/bin/env bash



echo "Hello from execute_for_crab!"
echo "Here is the Input File PSet :"
python -c "import PSet; print PSet.process.dumpPython()"

echo "There are $# arguments : $*"

echo "Now executing"


beforemsg=""
aftermsg=""
for i in "$@"
do
case $i in
    Before=*)
    beforemsg="${i#*=}"
    ;;
    After=*)
    aftermsg="${i#*=}"
    ;;
esac
done

echo "================= $beforemsg ===================="
python execute_for_crab.py
echo "================= $aftermsg ===================="
