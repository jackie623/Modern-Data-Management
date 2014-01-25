#!/bin/sh

echo "My current directory is:"
pwd

echo "Here are the files in my current directory:"
ls -la

echo "We're now changing to the sample-files directory"
cd sample-files

echo "These are my sample files:"
ls -la

echo "What's in this file?"
cat test-csv.csv

echo "For large files, you can use 'head' or 'less' to read the file"
echo "How does the 'head' command work?"
man head

head -3 test-csv.csv
less test-csv.csv # Press Q to exit less window

echo "Now that we're comfortable looking at the content of files, let's try moving files around"

cp test-csv.csv large-csv-file.csv
ls -la

mv large-csv-file.csv plant-experiment-data.csv
ls -la 

echo "Keep in mind that the mv command moves the file to the new path AND overwrites whatever's there. Unlike doing this in the finder window, this operation can't be undone, so BE CAREFUL with this one."

echo "If you're careful with the mv command, be extra extra careful with the rm command. It doesn't move a file to the trash, instead the file gets removed instantly and PERMENANTLY. Use the -r option do delete an entire directory and its contents." 
rm test-csv.csv

echo "Up till now, we've covered tasks that can be done, more or less from the Finder window (or he File Explorer, if you're on windows). Let's try something really powerful that you can't easily do from a built-in UI tool."

echo "Let's search our current directory for files with names matching a specific criteria. We'll start by finding all the csv files in our current directory. Note that the first argument we pass to the find command is a single dot (.). This dot is shorthand for the current directory."
find . -name "*.csv" # control-C kills a command that's taking too long

echo "Just to prove that this operation will search an entire directory (including subdirectories), let's move up a directory and see if we can still find the csv file. Note that the double dot (..) is shorthand for the directory above the current directory."
cd ..
find . -name "*.csv"

echo "Great! Now let's try to search on files that contain a specific word inside them. The options -ir that I've passed in to grep do two things - the i allows grep to match the text regardless of case, and the r tells grep to search the file we've provided - current directory (.) in this case - and any files below it in the directory tree. "
grep -ir "clementine" .

echo "Find all the csv files in the current directory only. To do this, we'll introduce a really powerful tool called the pipe (|). " 
ls -la | grep "\.csv"

echo "For good measure, let's leave this directory in the state that we found it, by renaming our csv file back to test-csv."

mv sample-files/plant-experiment-data.csv  sample-files/test-csv.csv


echo "Exercise 1: I'm a scientist doing some groundbreaking work on plant growth. Because you asked, and because I'm required to by law, I'm sharing my data with you. You've already downloaded my tarballed data set into your current working directory - it's called jackies-plant-data.tar. Uncompress the data, take a look at what's in there, and get a sense of what I did in my experiment."

echo "Exercise 2: You've noticed that my data isn't as organized as you like. In particular, you'd like the csvs to be split out into their own folders by date. Fix this problem using nothing but the shell." 