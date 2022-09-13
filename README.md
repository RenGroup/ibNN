# ibNN
Interpretable bionic neural network

General introduction</br>
ibNN is a novel neural network which simulates the structure of human signaling and gene regulatory network, incorporates existing biological knowledge, and learns the molecular relations from single cell RNA-seq data. The core of the network is built upon the conversion of the adjacency matrix of the topological directed graph of signaling network and TF-target relations to the initial weight matrices of ibNN. Therefore, each node in the network has been assigned an explicit name of the genes, and the trained weight matrices can be converted back to the relations between molecules, which provides clear intepretation of the meaning of the trained network.

The core design of ibNN</br>
The core of ibNN is the conversion of the directed graph of biological networks to the initial weight matrices of neural network. The original idea of the conversion is as the following figure:
<img width="1273" alt="adjacencyMatrixConvertion" src="https://user-images.githubusercontent.com/109563761/189895896-2aee0246-b5b4-49f5-99da-e72b0e2a000a.png">

Figure 1. The conversion of graphical knowledges of regulation to weight matrices of neural network. (A) A simple directed graph of three proteins. Protein A activates protein B, and protein B activates protein C. This type of graph is often seen in the signaling cascade, e.g. the PPrel of KEGG database. (B) The adjacency matrix of A. (C) A two-layer directed graph rearranged from A. (D) The adjacency matrix of the directed graph C with three vertices in a general format. (E) The layout of the input and middle layers of ibNN. For demonstration purpose, only three genes: A, B and C were shown. (F) The formula to calculate the input signal to the middle layer. The Winput_middle was derived from the transpose of the adjacency matrix in (D).

