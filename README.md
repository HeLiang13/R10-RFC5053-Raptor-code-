# R10-RFC5053-Raptor-code-
Provide encoder and decoder of R10
This projects provides encoder and decoder of Raptor codes, which is detailed in the standard file "RFC5053".
It is written with MATLAB.

Four decoding methods are included. 
1. The stadard decoding algorithm in RFC5053;
2. Random choose for pivot row in decoding Phase 1.
3. Kim's method for choosing pivot row in Phase 1, which can be found in paper of DOI 10.1109/LCOMM.2008.080599;
4. Zhang's algorithm which has pre-processed the decoding matrix, which is in paper of DOI 10.1109/ICMTMA.2010.531.

In Kim's and Zhang's , they both declare that their algorithm is better than the standard one in terms of decoding time. But in my simulation, Kim's is the same as the random method, and Zhang's is worse than the standard one. Looking forward to your correction.
My simulation results are in file "time_u_overhead_1.01.mat" and "time_u_overhead_1.05.mat"(here overhead is literally the decoding overhead),  corresponding to Kim's method，including 1.,2., and 3. method shown above. And file "time_u_Zhang_overhead_1.05.mat" is corresponding to Zhang's method.

Run "encode.m' and "decode.m", you get the results of encoder where you can set the block length "K", and then you get the results of decoder. If you want to run simulation for decoding time, use "Main_simulation.m" and uncomment relevant line in "encode.m' and "decode.m" in order to modify the two ".m" file into function file.

-----------------------------------------------------------------------------------------
本项目提供RFC 5053中Raptor码的编译码程序。
译码算法包括：
1.标准算法
2.随机选取行的方法
3.Kim的方法
4.Zhang的方法。

3.见文献An Efficient Algorithm for ML Decoding of Raptor Codes over the Binary Erasure Channel。
4.见文献An Improved Algorithm of 3GPP MBMS Raptor codes。
1.2.3.4的方法集合在decode.m中。

time_u_overhead_1.01.mat 和 time_u_overhead_1.05.mat 
给出的是1-3的结果，其overhead如文件名，K值设定基于3对应的文献，即K_list = [1000:1000:7000, 8192];
time_u_Zhang_overhead_1.05.mat给出的4的结果，K值设定基于4对应的文献，即K_list = 200:200:2000;

"encode.m' 和 "decode.m" 是编译器和译码器，可以先运行前者，再运行后者，其中可以设置信息源块的长度K，为了仿真方便，信息符号的大小我设置为1bit，可以自行改变。若要进行译码时间的蒙特卡洛仿真，则运行 "Main_simulation.m"， 但同时取消"encode.m' 和 "decode.m"中的相关注释，使这两个文件变为函数形式以便调用。

