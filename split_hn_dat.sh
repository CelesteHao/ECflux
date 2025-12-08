#!/bin/bash
shopt -s nocaseglob    # <-- enable case-insensitive filename matching

# BASEDIR="$HOME/Split bash/"
OUTDIR="/mnt/c/split/hn_1day/"  # Save outputs in WSL home

mkdir -p "$OUTDIR"

# Loop through each directory listed in dir_list.txt
while read -r DATDIR; do
    echo "Processing $DATDIR ..."

    # Make sure the directory exists
    if [ ! -d "$DATDIR" ]; then
        echo "Directory $DATDIR does not exist, skipping."
        continue
    fi

	for FILE in "$DATDIR"/*_TS10Hz.dat; do
		echo "Processing $FILE"

		head -n 4 "$FILE" > header.tmp

		awk -F, -v outdir="$OUTDIR" '
		NR > 4 {
			split($1, dt, " ");
			date = substr(dt[1], 2);
			out = outdir "/Hn1_" date "_TS.dat";

			if (!(out in seen)) {
				if (!system("[ -f " out " ]")) {
					seen[out] = 1;
				} else {
					system("cat header.tmp > " out);
					seen[out] = 1;
				}
			}

			print $0 >> out;
		}' "$FILE"

		rm -f header.tmp
	done

	for FILE in "$DATDIR"/*_Met1min.dat; do
		echo "Processing $FILE"

		head -n 4 "$FILE" > header.tmp

		awk -F, -v outdir="$OUTDIR" '
		NR > 4 {
			split($1, dt, " ");
			date = substr(dt[1], 2);
			out = outdir "/Hn1_" date "_MET.dat";

			if (!(out in seen)) {
				if (!system("[ -f " out " ]")) {
					seen[out] = 1;
				} else {
					system("cat header.tmp > " out);
					seen[out] = 1;
				}
			}

			print $0 >> out;
		}' "$FILE"

		rm -f header.tmp
	done

	for FILE in "$DATDIR"/*_Soil_1min.dat; do
		echo "Processing $FILE"

		head -n 4 "$FILE" > header.tmp

		awk -F, -v outdir="$OUTDIR" '
		NR > 4 {
			split($1, dt, " ");
			date = substr(dt[1], 2);
			out = outdir "/Hn1_" date "_SOIL.dat";

			if (!(out in seen)) {
				if (!system("[ -f " out " ]")) {
					seen[out] = 1;
				} else {
					system("cat header.tmp > " out);
					seen[out] = 1;
				}
			}

			print $0 >> out;
		}' "$FILE"

		rm -f header.tmp
	done

done < dir_list.txt   # <- reads the file line by line
echo "Done processing all .dat files."
