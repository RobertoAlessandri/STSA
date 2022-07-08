# Short-Time Spectral Attenuator and process reversion for musical effects

## Description

In this work, we implemented a STSA and tested its main configurations: as adapting Wiener Filter, as Spectral Attenuator and as Power Spectral Attenuator. 
After that we tried to see if it was possible to get a musical effect out of it. And we did it by inverting the suppresion rule and by applying the algorithm to clean signals instead of those with global degradations. In the OLA phase, instead of reconstructing the signal, we resynthesized it with the use of sinusoidal oscillators.


## Prerequisities

* To have installed on your local machine, the libraries imported at the beginnig of the notebooks if you don't use Colab
* MatLab 

## How to use

### Python Case

* Upload in your Colab Notebook all the [Audio Files](https://github.com/RobertoAlessandri/STSA/tree/main/Audio%20Files);

* Save a copy of both notebooks in your Drive

* Use Hiss folder if you want to test the [STSA](https://github.com/RobertoAlessandri/STSA/blob/main/STSA.ipynb). Try to experiment by changing the values of the "a, b, c" parameters

* Use Clean folder if you want to test the [Music Effect](https://github.com/RobertoAlessandri/STSA/blob/main/Reverse_STSA_MusicFX.ipynb). Try to experiment by changing the values of the "a, b, c" parameters.

### MatLab Case

* Download All the MatLab and Audio files. Place them all in the same folder without any sub-folder division

* Use audio files that start with "Source" if you want to test the [STSA](https://github.com/RobertoAlessandri/STSA/blob/main/MatLab%20Version/STSA.m). Try to experiment by changing the values of the "a, b, c" parameters

* Use audio files that end with "clean"  or "dehissed "if you want to test the [Music Effect](https://github.com/RobertoAlessandri/STSA/blob/main/MatLab%20Version/Reverse_STSA_MusicFX.m). Try to experiment by changing the values of the "a, b, c" parameters.

## Source Codes

All the codes are based on the theory explained by Godsin and Reiner in their papers and books of late 90's. TO DO PUTE REFERENCES AND THEN HYPERLINKS
(Also the audio files come from their site). TODO HYPERLINKS

Here I'll talk about Python files, but the structure for the MatLab ones is the same.
The only difference is that in MatLab functions are defined in different files. But that won't change the execution.
The scripts has the following structure:

* Imports of requried libraries;

* Methods definitions. The Music Effect file in this point has also the definition of a sinusoidal oscillator and of the resynthesizer.

* Loading of the files

* Initialize the parameters for the Short-Time Fourier Transform (STFT)

* Noise Estimation

* Suppression Rule. Here in the Music Effect file we resynthesized the signal

* Train dataset generation -> To generate the dataset for the training, first we randomly select 300 sentences from the TIMIT/TRAIN dataset. And since we want to reach a dataset size of 6000 datapoints, we locate 1 sentence in 30 different DOAs;

* Final Results and Plots.

## Presentation

A detailed description of our work can be find here: **TODO** put hypertext to the slides (?)

## Contacts

* Alessandri Roberto: roberto.alessandri@mail.polimi.it
* Cicognani Roberto Leone: robertoleone.cicognani@mail.polimi.it
