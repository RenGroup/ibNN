import numpy
import scipy.special
import re
import random
import sys
import os
import matplotlib.pyplot as plt
import time
import subprocess #for counting the lines of input .csv
import argparse
import pandas as pd #for batch output

# create parser object
parser = argparse.ArgumentParser(description = "Train ibNN and impute gene expressions")

# defining arguments for parser object
parser.add_argument("-d", "--dir", type = str, nargs = 1,
                    metavar = "path_to_csv", default = None,
                    help = "The path to the input csv.")
parser.add_argument("-i", "--filename", type = str, nargs = 1,
                    metavar = "input_file_name", default = None,
                    help = "The name of the input csv file.")
parser.add_argument("-e", "--nRounds", type = int, nargs = 1, #or number of epochs
                    metavar = "Max_num_training", default = 20,
                    help = "Training will end if flag turns to 3 or rounds of training reach this number")
parser.add_argument("-l", "--learningRate", type = float, nargs = 1,
                    metavar = "Learning_rate", default = 0.02, #set by experiences. Tried 0.1, 0.01, 0.005, and 0.001
                    help = "set by experiences.")
parser.add_argument("-c", "--nTrain", type = int, nargs = 1,
                    metavar = "Num_cell_for_training", default = -1,
                    help = "Num of cells for training, default '-1' for auto detect (2/3 of total cell number, may take time to count)") #this greatly effects runtime
parser.add_argument("-b", "--batchSize", type = int, nargs = 1,
                    metavar = "Num_cell_in_one_batch", default = 1,
                    help = "The mean of the errors of the cells in one batch will be used for back propagation")

# parse the arguments from standard input
args = parser.parse_args()

if args.dir != None:
    folder=str(args.dir[0])
    if(os.path.exists(folder)):
        print("The input dir exists.")
    else:
        exit("Dir ",folder," does not exist.")
    if(not folder.endswith('/')):
        folder=folder+"/"
else:
    print("Dir to the input csv file is required.")

if args.filename != None:
    #folder=str(args.dir[0])
    filename=str(args.filename[0])
    if(os.path.exists(folder+filename)):
        print("The input file exists.")
    else:
        if(not filename.endswith('.csv')) & (not filename.endswith('.txt')):
            exit("Please check the file format, or the suffix (need to be .csv or .txt)")
        else:
            exit("File ",folder+filename," does not exist.")
else:
    exit("The input csv or txt file is required.")

#Create the log file
folder_log = re.sub("masked/","",folder)+"log/"
if(not os.path.isdir(folder_log)):
    print(folder_log," does not exist, creating the folder for log files.")
    os.mkdir(folder_log)
file_log = "log_"+filename+".txt"
out_log = open(folder_log+file_log,'w')
print("The log file can be found at",folder_log+file_log)

if (type(args.nRounds) != int): #if use stdin value, then args.nRounds is Namespace; if use default value, then it's int
    #print(type(args.nRounds))
    n_rounds=args.nRounds[0]
else:
    n_rounds=args.nRounds

if (type(args.learningRate) != float):
    learning_rate = args.learningRate[0]
else:
    learning_rate = args.learningRate

if (type(args.nTrain) != int):
    n_train=args.nTrain[0]
else:
    n_train=args.nTrain

if (type(args.batchSize) != int):
    n_batch=args.batchSize[0]
else:
    n_batch=args.batchSize

if(n_train == 0):
    out_log.writelines("Error: entered zero cells for training. Please specify a positive number or use '-1' for auto detection.\n")
    exit("Error: entered zero cells for training. Please specify a positive number or use '-1' for auto detection.")
elif(n_train == -1):
    out_log.writelines("Selected auto detection of the cell numbers for training.\nStart counting the number of cells in input.\n")
    print("Selected auto detection of the cell numbers for training.")
    print("Start counting the number of cells in input.")
    cmd = "wc -l "+folder+filename
    cmdOut = subprocess.check_output(cmd,shell=True).decode("utf-8")
    cmdRegex = re.compile(r'\d+ ')
    lineCount = int(cmdRegex.search(cmdOut).group())
    n_train=int((lineCount-1)/3*2)
    out_log.writelines("There are "+str(lineCount-1)+" cells from the input, the first "+str(n_train)+" cells will be used for training.\n")
    print("There are "+str(lineCount-1)+" cells from the input, the first "+str(n_train)+" cells will be used for training.")

#if the num of cells for training is too low (<100), the possiblity of failing the mse filter increases.
#with more rounds (epochs) of training, the model may acquire sufficient accuracy.
if(n_train < 100) & (n_rounds < 40):
    out_log.writelines("Auto-adjusting triggered: training cell number is low ("+str(n_train)+"), n_rounds may not be sufficient ("+str(n_rounds)+").\nAdjusted n_rounds to be 40.\n")
    print("Auto-adjusting triggered: training cell number is low (",n_train,"), n_rounds may not be sufficient (",n_rounds,").")
    print("Adjusted n_rounds to be 40.")
    n_rounds = 40

#if n_train is too big, the run time will greatly affect the efficiency of ibNN since the accuracy may not further increase
#to accelerate the training of ibNN when n_train > 10000, batchSize will be automatically set to 1/1000 of n_train

if(n_batch == 1) & (n_train > 5000):
    n_batch = int(n_train/1000)
    out_log.writelines("Auto-adjusting triggered: training cell number is big ("+str(n_train)+"), batchSize will be adjusted to "+str(n_batch)+".\n")
    print("Auto-adjusting triggered: training cell number is big (",n_train,"), batchSize will be adjusted to",n_batch,".")

#adjust the chunk size and coefficient:
if(n_batch == 1):
    n_chunk = 1000
elif(n_batch > 1):
    c_chunk = int(1000/n_batch) # c_chunk is the number of times of n_batch within one chunk. n_chunk is set to close to 1000, always smaller
    n_chunk = c_chunk * n_batch # when 1000 % n_batch != 0, this number is smaller than 1000

out_log.writelines("The parameters are:\nn_rounds: "+str(n_rounds)+" learning_rate: "+str(learning_rate)+" n_train: "+str(n_train)+" n_batch: "+str(n_batch)+" n_chunk: "+str(n_chunk)+".\n")
print("The parameters are:\nn_rounds:",n_rounds,"learning_rate:",learning_rate,"n_train:",n_train,"n_batch:",n_batch,"n_chunk:",n_chunk)
#This is the global nn version. Modified the error calculation method: only the genes with >0 values in target are used

#neural network class definition
class neuralNetwork :
    #initialise the neural network
    def __init__( self, inputnodes, hiddennodes, outputnodes, learningrate, wMa_wih, wMa_who) :
        # set number of nodes in each input, hidden, output layer
        self.inodes = inputnodes
        self.hnodes = hiddennodes
        self.onodes = outputnodes

        #link weight matrices, wih and who
        self.wih = wMa_wih.T # the shape of the wih is (hidden,input)
        self.who = wMa_who.T # the shape of the who is (output,hidden)

        #learning rate
        self.lr = learningrate

        #activation function is the sigmoid function
        self.activation_function = lambda x: scipy.special.expit(x) #signmoid

        pass
    #the method to train the neural network
    def train(self, inputs_list, targets_list) :
        #convert inputs list to 2d array
        inputs = numpy.array(inputs_list, ndmin = 2).T
        targets = numpy.array(targets_list, ndmin = 2).T

        #calculate signals into hidden layer
        hidden_inputs = numpy.dot(self.wih, inputs)
        #calculate the signals emerging from hidden layer
        hidden_outputs = self.activation_function(hidden_inputs)

        #calculate signals into final output layer
        final_inputs = numpy.dot(self.who, hidden_outputs)
        #calculate the signals emerging from final output layer
        final_outputs = self.activation_function(final_inputs)
        #!IMPORTANT: setting zeros for the genes with 0 values in targets
        indx_pos = [k for k in range(len(targets)) if float(targets[k]) > 0]
        indx_zero = list(set(range(len(targets))) - set(indx_pos))
        final_outputs[indx_zero] = 0 #setting the genes with 0 value in targets to 0

        #output layer error is (target - actual)
        output_errors = targets - final_outputs
        self.output_errors = targets - final_outputs
        #hidden layer error is the output_errors, split by weights, recombined at hidden nodes
        hidden_errors = numpy.dot(self.who.T, output_errors)

        #update the weights for the links between the hidden and output layers
        self.who += self.lr * numpy.dot((output_errors * final_outputs * (1.0 - final_outputs)), numpy.transpose(hidden_outputs))

        #update the weights for the links between the input and hidden layers
        self.wih += self.lr * numpy.dot((hidden_errors * hidden_outputs * (1.0 - hidden_outputs)), numpy.transpose(inputs))

        pass
    #the method to train the neural network for batch input, taking the median of the batch errors. "inputs" is a 2-d array
    def train_batch(self, inputs_list, targets_list) :
        #convert inputs list to 2d array
        inputs = numpy.array(inputs_list, ndmin = 2).T
        targets = numpy.array(targets_list, ndmin = 2).T

        #calculate signals into hidden layer
        hidden_inputs = numpy.dot(self.wih, inputs)
        #calculate the signals emerging from hidden layer
        hidden_outputs = self.activation_function(hidden_inputs)

        #calculate signals into final output layer
        final_inputs = numpy.dot(self.who, hidden_outputs)
        #calculate the signals emerging from final output layer
        final_outputs = self.activation_function(final_inputs)

        #!IMPORTANT: setting zeros for the genes in outputs with 0 values in targets. Loop for cells
        for j in range(targets.shape[1]): #each column is a cell now
            indx_pos = [k for k in range(targets.shape[0]) if float(targets[k,j]) > 0]
            indx_zero = list(set(range(targets.shape[0])) - set(indx_pos))
            final_outputs[indx_zero,j] = 0 #setting the genes with 0 value in targets to 0

        #output layer error is (target - actual). Note this is a 2-d array now; each columan is a cell, each row is a gene
        output_errors = numpy.array(numpy.mean(targets - final_outputs,axis = 1), ndmin = 2).T #output_errors is a 1-d array
        self.output_errors = numpy.array(numpy.mean(targets - final_outputs,axis = 1), ndmin = 2).T
        #hidden layer error is the output_errors, split by weights, recombined at hidden nodes
        hidden_errors = numpy.dot(self.who.T, output_errors)

        #update the weights for the links between the hidden and output layers
        self.who += self.lr * numpy.dot((output_errors * final_outputs * (1.0 - final_outputs)), numpy.transpose(hidden_outputs))

        #update the weights for the links between the input and hidden layers
        self.wih += self.lr * numpy.dot((hidden_errors * hidden_outputs * (1.0 - hidden_outputs)), numpy.transpose(inputs))

        pass
    #query the neural network
    def query(self, inputs_list) :
        #convert input list to 2d array
        inputs = numpy.array(inputs_list, ndmin = 2).T

        #calculate signals into hidden layer
        hidden_inputs = numpy.dot(self.wih, inputs)
        #calculate the signals emerging from hidden layer
        hidden_outputs = self.activation_function(hidden_inputs)

        #calculate signals into final output layer
        final_inputs = numpy.dot(self.who, hidden_outputs)
        #calculate the signals emerging from final output layer
        final_outputs = self.activation_function(final_inputs)
        return final_outputs
        pass

#load weight matrix
wMa_pprel_file = open ("/Users/chix/neuralNetwork/wMa.pprel.txt",'r')
wMa_pprel_list = wMa_pprel_file.readlines()
wMa_pprel_file.close()
wMa_pprel_col_geneID = wMa_pprel_list[0]
wMa_pprel_col_geneID = wMa_pprel_col_geneID[:-1].split("\t")[1:] #remove "\n" and split by "\t"
for i in range(1,len(wMa_pprel_list)):
    record = wMa_pprel_list[i][:-1]
    wMa_pprel_list[i] = numpy.asfarray(record.split('\t'))
wMa_pprel_all = numpy.asarray(wMa_pprel_list[1:])
wMa_pprel_row_geneID = wMa_pprel_all[:,0].astype(int) #store the row gene ID, convert the float type to int
wMa_pprel_all = wMa_pprel_all[:,1:]

wMa_tf_file = open ("/Users/chix/neuralNetwork/wMa.tf.txt",'r')
wMa_tf_list = wMa_tf_file.readlines()
wMa_tf_file.close()
wMa_tf_col_geneID = wMa_tf_list[0]
wMa_tf_col_geneID = wMa_tf_col_geneID[:-1].split("\t")[1:] #remove "\n" and split by "\t"
for i in range(1,len(wMa_tf_list)):
    record = wMa_tf_list[i][:-1]
    wMa_tf_list[i] = numpy.asfarray(record.split('\t'))
wMa_tf_all = numpy.asarray(wMa_tf_list[1:])
wMa_tf_row_geneID = wMa_tf_all[:,0].astype(int) #store the row gene ID, convert the float type to int
wMa_tf_all = wMa_tf_all[:,1:]

#since all the cells have the same order of gene ids, the indices of the input can be processed before the loops
sc_data_file = open (folder+filename,'r')
sc_data_line = sc_data_file.readline() #read the first line, the gene ids
sc_data_line = re.sub("[\n\r]","",sc_data_line) #usage: re.sub(pattern, replacement, string)
geneID = sc_data_line.split(',')[1:] #removed the first "" element
sc_data_file.close()

#Note: there might be duplicates in geneID, which causes differences between local network version and global network version

#intersect input and weight matrix
##for inputs and wih
xy, xy_ind1, xy_ind2 = numpy.intersect1d(geneID, wMa_pprel_row_geneID, return_indices = True)
indx_input = xy_ind1 #for every cell, this is the index of the input

#the weight matrices still need to be subsetted, but only for one time
wMa_pprel_sub = wMa_pprel_all[xy_ind2,] #!the wih of the global nn
wMa_pprel_row_geneID_sub = wMa_pprel_row_geneID[xy_ind2]#the gene id of the wih of the global nn

##for who and targets
xy, xy_ind1, xy_ind2 = numpy.intersect1d(geneID, wMa_tf_col_geneID, return_indices = True)
indx_target = xy_ind1 #the indices of the target genes. it's the same for every cell. The duplicates will be dropped, only the first ocurrence of match will be retained
wMa_tf_sub = wMa_tf_all[:,xy_ind2] #!the who of the global nn
wMa_tf_col_geneID_sub = numpy.asarray(wMa_tf_col_geneID)[xy_ind2] #the gene id of the who of the global nn

#use global network rather than local network to avoid creating and updating the network for every cell
#number of input, hidden and output nodes
input_nodes = wMa_pprel_sub.shape[0] #although still a subset of the global network, it doesn't need to created for every single cell
hidden_nodes = wMa_pprel_sub.shape[1]
output_nodes = wMa_tf_sub.shape[1]

#create instance of neural network
n = neuralNetwork(input_nodes,hidden_nodes,output_nodes,learning_rate, wMa_pprel_sub, wMa_tf_sub)

#initialize array for internal cell-level MSE
updated_genes = numpy.array([])
arr_squared_all = numpy.array([])

#n_rounds = 20 #max number of epochs of rounds. Usually less than 20 rounds are needed
flag = 1; #1 to continue training; 2 to stop
out_log.writelines("Finished initialization. Start training the neural network.\n")
print("Finished initialization. Start training the neural network.")

for j in range(n_rounds): #for each round, the file will be opened once.
    #sc_data_file = open (folder+filename,'r')
    #sc_data_line = sc_data_file.readline() #read the first line, the gene ids. Nothing to do, move on
    n_cells=0
    inputs = numpy.array([])
    arr_squared=numpy.array([])
    start = time.time()
    r_chunk = 0 #round of chunks
    n_debug = 0
    for chunk in pd.read_csv(folder+filename,chunksize=n_chunk): #the header is automatically ignored
        chunk=chunk.to_numpy()
        #print(chunk.shape)
        cell_id=chunk[:,0]
        chunk=chunk[:,1:] #keep the expr values only, rm the cell names
        r_chunk = r_chunk + 1
        if(n_batch == 1):
            indx_cell = -1
            #train when n_batch ==1
            while (n_cells < n_train) & (indx_cell+1 < chunk.shape[0]):
                indx_cell = n_cells % n_chunk #n_cells start from zero
                n_cells = n_cells +1
                inputs = numpy.array(chunk[indx_cell,:])[indx_input] #input genes, including zeros. Less than 3500. n_cells start from zero
                targets = numpy.array(chunk[indx_cell,:])[indx_target] #!IMPORTANT: since indx_target is for all cells, there are genes whose expr are zeros in this array
                tmp_geneid = numpy.array(geneID)[indx_target] #the ids matching targets
                indx_pos = [k for k in range(len(targets)) if float(targets[k]) > 0] # the genes with positive values contribute to the updates of the matrices
                tmp_geneid = tmp_geneid[indx_pos] #select the positive genes from "targets"
                updated_genes = numpy.union1d(updated_genes, tmp_geneid) #this is useful when output the imputation results. Only trained genes will be output

                inputs = numpy.asfarray(inputs)
                targets = numpy.asfarray(targets)
                inputs = numpy.log2(inputs +1) #log transformation
                targets = numpy.log2(targets +1) #log transformation
                scaleFac= numpy.max(targets) # targets includes all the genes with positive values. Scale the data to [0,1].
                targets = targets/scaleFac
                inputs = inputs/scaleFac

                #train the neural network
                epochs = 1 #better MSE, but not much. Need to leverage with time consumption. Set to 1 if the number of training cell > 1000
                for e in range(epochs):
                    n.train(inputs,targets)

            #calculate mse
            while(n_cells >= n_train) & (indx_cell+1 < chunk.shape[0]) & (n_cells < lineCount-1): #lineCount - 1 is the total number of cells in the input
                indx_cell = n_cells % n_chunk
                n_cells = n_cells +1
                inputs = numpy.array(chunk[indx_cell,:])[indx_input] #input genes, including zeros. Less than 3500
                inputs = numpy.asfarray(inputs)
                inputs = numpy.log2(inputs +1) #log2 transformation
                targets = numpy.asfarray(chunk[indx_cell,:])[indx_target]
                targets = numpy.log2(targets + 1)
                scaleFac= numpy.max(targets)
                targets = targets/scaleFac
                targets = numpy.array(targets, ndmin = 2).T
                inputs = inputs/scaleFac

                #query the neural network
                output_query = n.query(inputs)

                 #!IMPORTANT: setting zeros for the genes with 0 values in targets
                indx_pos = [k for k in range(len(targets)) if float(targets[k]) > 0]
                indx_zero = list(set(range(len(targets))) - set(indx_pos))
                output_query[indx_zero] = 0 #setting the genes with 0 value in targets to 0

                df_subtr = (output_query - targets)*scaleFac #need to scale back
                gene_num = len(indx_pos) #this is the number of positive genes
                df_squared = numpy.sum(numpy.power(df_subtr,2))/gene_num
                #print(df_squared)
                arr_squared = numpy.append(arr_squared,df_squared)

        elif(n_batch >1):
            indx_cell = -1
            #train when n_batch >1
            #print(n_cells,n_train,n_chunk,chunk.shape[0])
            while (n_cells < n_train) & (n_cells < (r_chunk-1)*n_chunk+chunk.shape[0]):
                indx_cell = n_cells % n_chunk #the index of the start of the batch
                indx_end = indx_cell + n_batch  # the index of the end of the batch
                if(n_cells + n_batch -1 >= n_train): #the end of the batch has to be adjusted when reaches the n_train or chunk.shape[0]
                    indx_end = n_train % n_chunk
                    n_cells = n_train
                elif(indx_end -1 > chunk.shape[0]):
                    indx_end = chunk.shape[0]+1
                    n_cells = n_cells + chunk.shape[0] - indx_cell + 1
                else:
                    n_cells = n_cells + n_batch
                #print(indx_cell,indx_end)
                #print("Position 1:",numpy.mean(indx_input))
                inputs = numpy.asfarray(chunk[indx_cell:indx_end,indx_input])
                targets = numpy.asfarray(chunk[indx_cell:indx_end,indx_target]) #this is only for the calculation of scaleFac
                #print(inputs.shape,targets.shape)
                tmp_geneid = numpy.array(geneID)[indx_target] #the ids matching targets
                col_max = numpy.max(targets,axis = 0) #the max number of counts of each gene, to find non-zero genes
                indx_pos = [k for k in range(len(col_max)) if float(col_max[k]) > 0] # the genes with positive values contribute to the updates of the matrices
                tmp_geneid = tmp_geneid[indx_pos] #select the positive genes
                updated_genes = numpy.union1d(updated_genes, tmp_geneid) #this is useful when output the imputation results. Only trained genes will be output

                inputs = numpy.log2(inputs +1) #log transformation
                targets = numpy.log2(targets +1) #log transformation
                scaleFac= numpy.max(targets,axis = 1) # targets includes all the genes with positive values. Scale the data to [0,1].
                #print(len(scaleFac),inputs.shape,targets.shape)
                #print(numpy.around(scaleFac,4),inputs.shape,targets.shape)
                #exit()
                for i in range(inputs.shape[0]): #scale the inputs by its own max expr values
                    targets[i,:] = targets[i,:]/scaleFac[i]
                    inputs[i,:] = inputs[i,:]/scaleFac[i]
                #print(numpy.mean(inputs[0,:]),numpy.mean(targets[0,:]))
                #exit()
                #train the neural network
                epochs_inner = 1 #>1 may get better MSE, but not much. Need to leverage with time consumption. Recommend "1" if the number of training cell > 1000
                for e in range(epochs_inner):
                    if(inputs.shape[0]> 1):
                        n.train_batch(inputs,targets)
                    elif(inputs.shape[0] == 1):
                        n.train(inputs,targets)
                    #n_debug = n_debug + 1
            #print(n_debug)
            #print(n_cells,n_train,n_chunk,chunk.shape[0])
            #calculate mse
            while(n_cells >= n_train) & (indx_cell+1 < chunk.shape[0]) & (n_cells < lineCount-1): #lineCount - 1 is the total number of cells in the input
                indx_cell = n_cells % n_chunk
                #print(indx_cell,n_cells,cell_id[indx_cell])
                n_cells = n_cells +1
                #print("Position 2:",numpy.mean(indx_input))
                inputs = numpy.array(chunk[indx_cell,:])[indx_input] #input genes, including zeros. Less than 3500
                inputs = numpy.asfarray(inputs)
                inputs = numpy.log2(inputs +1) #log2 transformation
                targets = numpy.asfarray(chunk[indx_cell,:])[indx_target]
                targets = numpy.log2(targets + 1)
                scaleFac= numpy.max(targets)
                targets = targets/scaleFac
                targets = numpy.array(targets, ndmin = 2).T
                inputs = inputs/scaleFac

                #print(numpy.mean(inputs),numpy.mean(targets))
                #exit()
                #print(numpy.around(scaleFac,4),chunk.shape,inputs.shape,targets.shape)
                #exit()

                #query the neural network
                output_query = n.query(inputs)
                #print(output_query[0:5])
                #exit()

                 #!IMPORTANT: setting zeros for the genes with 0 values in targets
                indx_pos = [k for k in range(len(targets)) if float(targets[k]) > 0]
                indx_zero = list(set(range(len(targets))) - set(indx_pos))
                output_query[indx_zero] = 0 #setting the genes with 0 value in targets to 0

                df_subtr = (output_query - targets)*scaleFac #need to scale back
                gene_num = len(indx_pos) #this is the number of positive genes
                df_squared = numpy.sum(numpy.power(df_subtr,2))/gene_num
                #print(df_squared)
                arr_squared = numpy.append(arr_squared,df_squared)
    #store the trained weight matrices in case the training stopped because of over-fitting
    if(j == 0):
        prev_r_2_wih = n.wih
        prev_r_2_who = n.who
    elif(j == 1):
        prev_r_1_wih = n.wih
        prev_r_1_who = n.who
    elif(j > 2):
        prev_r_2_wih = prev_r_1_wih
        prev_r_2_who = prev_r_1_who
        prev_r_1_wih = n.wih
        prev_r_1_who = n.who

    #store the MSE. If you want to plot the MSEs for epochs, output the arr_squared_all at the end of the script
    if(len(arr_squared_all) == 0):#if they are void, initialize
        arr_squared_all = arr_squared
        out_log.writelines("The median mse now is: "+str(numpy.around(numpy.median(arr_squared),5))+".\n")
        print("The median mse now is:",numpy.around(numpy.median(arr_squared),5),".")
    else:
        arr_squared_all = numpy.vstack((arr_squared_all,arr_squared))

    #print and check the changes of median MSE
    if(j > 0):
        median_this = numpy.median(arr_squared)
        median_prev = numpy.median(arr_squared_all[j-1,:])
        #print("Shape of mse stack:",arr_squared_all.shape)
        out_log.writelines("prevMSE: "+str(numpy.around(median_prev,5))+" thisMSE: "+str(numpy.around(median_this,5))+" changePerc: "+str(numpy.around(abs((median_this - median_prev)/median_prev),8))+" lr: "+str(learning_rate)+".\n")
        print("prevMSE: "+str(numpy.around(median_prev,5))+" thisMSE: "+str(numpy.around(median_this,5))+" changePerc: "+str(numpy.around(abs((median_this - median_prev)/median_prev),8))+" lr: "+str(learning_rate)+".")
        #exit()
        if((median_this < 0.4) & (abs(median_prev - median_this)/median_prev < 0.001)) or (median_this - median_prev > 0):
            if(flag == 1):
                flag = flag +1
            elif (flag == 2):
                flag = flag +1
                out_log.writelines("Reached flag == 3, ending the training process at round "+str(j+1)+".\n")
                print("Reached flag == 3, ending the training process at round ",j+1,".")
                out_log.writelines("The trained weight matrices have been rolled back by two rounds.")
                print("The trained weight matrices have been rolled back by two rounds.")
                #roll back the trained weight matrices
                n.wih = prev_r_2_wih
                n.who = prev_r_2_who
                break
        else:
            flag = 1
        out_log.writelines("Flag: "+str(flag)+"\n")
        print("Flag:",flag)


    end =  time.time()
    out_log.writelines("Training and internal MSE calculation for round "+str(j+1)+" takes "+str(numpy.around(end-start,2))+ "s.\n\n")
    print("Training and internal MSE calculation for round ",j+1," takes ", numpy.around(end-start,2), "s.\n")


if(median_this > 1000): #set to a big number (such as 10000) to ignore
    out_log.writelines("Warning: internal control of MSE did not pass the filter (median MSE:"+str(median_this)+").\n Imputation can be risky, exit ibNN.\n Modify filter at line 512 to ignore this.\n")
    exit("Warning: internal control of MSE did not pass the filter (median MSE:"+str(median_this)+").\n Imputation can be risky, exit ibNN.\n Modify filter at line 512 to ignore this.")

out_log.writelines("Done training, start imputation.\n")
print("Done training, start imputation.")
start = time.time()
#impute cell by cell
#open the output file
out_folder=re.sub("masked/","",folder)+"imputed/"
#check if the required path already exists or not. If not, create one.
if(not os.path.isdir(out_folder)):
    out_log.writelines(out_folder+" does not exist, creating one.\n")
    print(out_folder," does not exist, creating one.")
    os.mkdir(out_folder)
out_file = "imputed_"+filename
out_imputed = open(out_folder+out_file,'w')

out_log.writelines("The output imputed expr file can be found at "+out_folder+out_file+"\n")
print("The output imputed expr file can be found at",out_folder+out_file)
#open the original data file
sc_data_file = open (folder+filename,'r')
sc_data_line = sc_data_file.readline() #read the first line, the gene ids. This is not the header of the output file

##for imputation, the output genes are from updated_genes who has been trained. This has been fixed, no need to do this for each cell
xy, xy_ind1, xy_ind2 = numpy.intersect1d(updated_genes, wMa_tf_col_geneID_sub, return_indices = True)
wMa_tf_updated = n.who.T[:,xy_ind2] #call for the trained who from "n"
wMa_tf_col_geneID_sub2 = wMa_tf_col_geneID_sub[xy_ind2] #subset of wMa_tf_col_geneID_sub to use only updated_genes

#use global network rather than local network to avoid creating and updating the network for every cell
#number of input, hidden and output nodes
input_nodes = wMa_pprel_sub.shape[0] #although still a subset of the global network, it doesn't need to create for every single cell
hidden_nodes = wMa_pprel_sub.shape[1]
output_nodes = wMa_tf_updated.shape[1]

#create instance of neural network
n2 = neuralNetwork(input_nodes,hidden_nodes,output_nodes,learning_rate, n.wih.T, wMa_tf_updated)
#the gene ids do not change for wih, therefore directly call the wih from "n"
#the output weight matrix has been changed according to "updated_genes"

header = ",".join(map(str,wMa_tf_col_geneID_sub2))
header = ","+header+"\n" #add a void element

#output the header
out_imputed.writelines(header)
out_imputed.close()
sc_data_file.close()

n_cell = 0
n_totalCell = 0
for chunk in pd.read_csv(folder+filename,chunksize=1000):
    chunk=chunk.to_numpy()
    cell_name = chunk[:,0]
    chunk=chunk[:,1:] #keep the expr values only, rm the cell names
    n_cell = chunk.shape[0] #the size of the chunk. May not equal to 1000 for the last loop
    n_totalCell = n_totalCell + n_cell
    #print(n_cell)
    inputs = numpy.asfarray(chunk)[:,indx_input]
    inputs = numpy.log2(inputs +1) #log2 transformation
    targets = numpy.asfarray(chunk)[:,indx_target] #this is only for the calculation of scaleFac
    targets = numpy.log2(targets + 1)
    scaleFac= numpy.max(targets,axis = 1) #An array of length n_cell
    for i in range(n_cell): #scale the inputs by the max of its own expr values
        inputs[i,:]=inputs[i,:]/scaleFac[i]
    #query the neural network
    output_query = n2.query(inputs).T# #only for dev, convert to n2 when paste into formal script
    for i in range(n_cell):
        output_query[i,:]=output_query[i,:]*scaleFac[i] #similar to the scaling of the inputs
    output_query = numpy.around(numpy.asfarray(output_query),4) #round numbers
    output=pd.DataFrame(output_query,index = cell_name) #convert to dataframe, add cell names
    output.to_csv(out_folder+out_file,mode='a',header=False)
    if(n_totalCell % 1000 == 0):
        out_log.writelines("Finished output "+str(n_totalCell)+" cells.\n")
        print("Finished output "+str(n_totalCell)+" cells.")

end=time.time()
out_log.writelines("Output imputation file takes "+str(end-start)+"s.\n")
print("Output imputation file takes ", end-start, "s")
#output weight matrices
#check if the required path already exists or not. If not, create one.
folder_wMa = re.sub("masked/","",folder)+"wma/"
print("Writing weight matrices.")
out_log.writelines("Writing weight matrices.\nWeight matrices can be found at "+folder_wMa+", begin with wih_ or who_\n")
print("Weight matrices can be found at ",folder_wMa,", begin with wih_ or who_")
out_file = re.sub(".csv",".txt",filename)
if(not os.path.isdir(folder_wMa)):
    out_log.writelines("Folder "+folder_wMa+" does not exist, creating one.\n")
    print("Folder "+folder_wMa+" does not exist, creating one.")
    os.mkdir(folder_wMa)
file = open(folder_wMa+"wih_"+out_file,"w")
header ="\t"+"\t".join(map(str,wMa_pprel_col_geneID))+"\n"
file.writelines(header)
wih = n.wih.T #the trained weight matrix of input to hidden (pprel)
for i in range(wih.shape[0]):
    joined_line = str(wMa_pprel_row_geneID_sub[i])+"\t"+"\t".join(map(str,wih[i,:]))
    joined_line = re.sub("[\[\]]","",joined_line)
    file.writelines(joined_line+"\n")
file.close()

file = open(folder_wMa+"who_"+out_file,"w")
header ="\t"+"\t".join(map(str,wMa_tf_col_geneID_sub))+"\n"
file.writelines(header)
who = n.who.T #the trained weight matrix of input to hidden (pprel)
for i in range(who.shape[0]):
    joined_line = str(wMa_tf_row_geneID[i])+"\t"+"\t".join(map(str,who[i,:]))
    joined_line = re.sub("[\[\]]","",joined_line)
    file.writelines(joined_line+"\n")
file.close()

out_log.writelines("The parameters are:\nn_rounds: "+str(j)+" learning_rate: "+str(learning_rate)+" n_train: "+str(n_train)+" n_batch: "+str(n_batch)+".\nDone output imputation result, ibNN finished.\n")
print("The parameters are:\nn_rounds:",j,"learning_rate:",learning_rate,"n_train:",n_train,"n_batch:",n_batch)
print("Done output imputation result, ibNN finished.")
out_log.close()
