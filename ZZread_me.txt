本文件夹提供RFC 5053中Raptor码的编译码程序。
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

