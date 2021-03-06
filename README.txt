Trina Rutz
April 7, 2018

Note: the program requires that the file path to the wave sample be inputted as a command line argument.

I began playing around with inputting the wave file in python pretty much as soon as it was announced. I had made some progress by the time the help files were posted, but I felt very unsure of why certain things needed to happen. Part of the problem is that I did not find any good resources explaining in detail how to read a wave file into python and what all the various steps meant, and so ‘findpeak.py’ was very useful to study and figure out exactly how the wave file was structured and how to read it in. Likewise, I also had mostly figured out the steps to solve the homework, but the ‘How-To’ was immensely helpful to confirm what I thought needed to happen and to use as a reference.

I found C code for the Goertzel filter from the University of Washington, but used other resources when in making my own goertzel function which resulted in one that didn’t produce correct values in every circumstance. Pat Rademacher helped me overhaul my goertzel function into one that worked, and my goertzel function is an altered form of his. It turned out that one of the resources I’d used when making my goertzel function was just entirely off. 

At first I did not normalize my results, and just called goertzel on every 160-element chunk and saved the results into an array. I figured out that one of two values was outputted from goertzel, and figured out that those numbers referred to which peak was higher. I changed the get_binary_from_peak function to save ’s’ or ‘m’ based on which result from the goertzel result was higher. I then mapped a 0 to each ’s’ and a 1 to each ‘m’. My results were weird until I normalized the wave sample values to be between -1 and 1, which helped the results stay more consistent.

Then I broke the array of binary values into 10-element chunks, and disregarded the 0th and 9th elements of every chunk after checking that the 0th element was a 0 and the 9th element was a 1. I then reversed the remaining chunk, converted it to an int, and converted that into a char through the binascii library.

My biggest problems while completing this assignment were just being very unsure of what to do until the ‘How-To’ and ‘findpeak.py’. After the program was completed I realized that it was very space inefficient, and would be improved by reading in the wave sample in 160-sample chunks, normalizing the values, calling goertzel on the entire chunk, and saving the correct character to a file before using the same array to store the next 160-sample chunk.
