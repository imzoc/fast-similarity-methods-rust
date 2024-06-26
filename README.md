# Pseudocode for src/scripts/generate_comparison_data.rs

### Constants
`K_range = [2, 3, 4]`  // this is subject to change.

`sequence_database_file_name = "sequences_1000.csv"` // this is my sequence database.

## Step 1: read metadata from the sequence database
The database file has three metadata lines at the top:
* “# base string length=1000”  // not always 1000.
* “# minimum edit distance=0  // not always 0.
* “# maximum edit distance=400  // not always 400.

__I'm having trouble parsing the numeric data (I should explain more).__

## Step 2: do computation on the database.
Each line in the database has 3 columns:
* Column 1 is a DNA sequence of variable length.
* Column 2 is a DNA sequence of variable length.
* Column 3 is a number indicating the true Levenshtein edit distance between the strings in columns 1 and 2.

We need to do computations for each line of data. Sp

## Step 3: do some statistics (not super relevant).
Step 4: save data to csv file
I’m saving average and variability data for each combination of “k” and “true Levenshtein edit distance”
The way I’m currently doing it:
rows = “true Levenshtein edit distance”
Columns = “kind of data (k={k})”
There will be 3 * |k_range| columns (not including a row header column).
There will be maximum edit distance - minimum edit distance rows (not including a column header row). 
