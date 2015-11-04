#!/bin/bash
./train 900 model_init.txt seq_model_01.txt model_01.txt
./train 900 model_init.txt seq_model_02.txt model_02.txt
./train 900 model_init.txt seq_model_03.txt model_03.txt
./train 900 model_init.txt seq_model_04.txt model_04.txt
./train 900 model_init.txt seq_model_05.txt model_05.txt
./test modellist.txt testing_data1.txt result1_900.txt >> acc_data1_900
./test modellist.txt testing_data2.txt result2_900.txt 
