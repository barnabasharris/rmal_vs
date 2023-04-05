module purge
module load r/recommended
module load gerun
export TZ="Europe/London"
export TMPDIR="/home/tcrnbgh/Scratch/tmp"

echo "TMPDIR=/home/tcrnbgh/Scratch/tmp" > /home/tcrnbgh/Scratch/visual_prominence/.Renviron
cp /home/tcrnbgh/Scratch/visual_prominence/.Renviron /home/tcrnbgh/.Renviron
