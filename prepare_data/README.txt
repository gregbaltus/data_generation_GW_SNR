Explanation of prepare_data.py
------------------------------

This function just copy the ../data_ratio3.0/noisy_GW/numpy into a folder call data_ratio3.0 
This function do atomatically that for each folder 
You have to specify the first folder, the last folder and the increment

python3 prepare_data.py --debut=3.0 --fin=16.0 --jump=0.5 --trainset=1

--trainset=1 if you want to have your data inside data_ratio3.0/trainset/ instead of just data_ratio3.0/ 
--testset=1 if you want to have your data inside data_ratio3.0/testset/ instead of just data_ratio3.0/


Note : this algo remove automatically the old folder call data_ratio3.0 for exemple
So if you don't have any folder it's not supprising you have this warning :
rm: cannot remove 'data_ratio4.5': No such file or directory
