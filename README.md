XPIRe: Private Information Retrieval for Everyone
=================================================

XPIRe allows a user to privately download an element from a database. This means that the database server knows that she has sent a database element to the user but does not know which one. The scientific term for the underlying protocol is Private Information Retrieval (PIR). This library is described and studied in the paper:

Carlos Aguilar-Melchor, Joris Barrier, Laurent Fousse, Marc-Olivier Killijian, "XPIRe: Private Information Retrieval for Everyone", available on the IACR eprints, http://eprint.iacr.org/2014/1025.

If you use our library, or a sub-part, such a NFLlib, please cite this paper on your work.

This project is closely related to another project available (soon?) at GitHub: NFLlib. The NTT-based Fast Lattice library, which allows fast cryptographic computations over ideal lattices. Adapting this project to a more recent version of NFLlib would provide a performance boost but would be a lot of work so we are not planning to do it immediately. 

*Important Note 1*: A PIR implementation for a specific application can be much simpler and more compact than this library. Indeed, much of the code of XPIRe is aimed to ensure that the library delivers good performance in a large span of applications *without user interaction* so that the user does not need to be an expert on PIR or cryptography to use it with good performance results. **If you want to use a PIR protocol for a very specific setting feel free to contact us for building up a collaboration !**

*Important Note 2*: For publication issues, a small part of the code is missing. From a technical point of view this correspond to the gaussian noise generator for LWE which is replaced by a uniform noise generator until some scientific results are published. Replacing our uniform noise generator with our gaussian noise generator does not impact performance in an observable way.

*Important Note 3*: This software cannot provide reliable privacy without more scrutiny on many details. We have tried to provide some resiliance to timing tests, but haven't tested them thoroughly. The random seed generation and pseudorandom generation use strong functions but we haven't done a thorough analysis of whether an idiotic fault is present in those critical sections of the code or not. No input verification is done by the server or client so many buffer overflows are potentially possible, etc. As is, the software *shows that privacy is possible* but cannot guarantee it against strong active adversaries (using timing attacks, introducing malformed entries, etc.) until it gets enough scrutiny.


Installation:
=============

Requirements: 
- 64-bits Linux OS: g++>=4.8, gcc>=4.8
- Apple OSX: [Xcode](https://itunes.apple.com/fr/app/xcode/id497799835?mt=12), [Xcode Command Line Tools](https://developer.apple.com/library/ios/technotes/tn2339/_index.html) and [MacPorts](https://www.macports.org/install.php).

Get a copy of the project with:
- git clone git@github.com:XPIRe-team/XPIRe.git
- or by downloading from https://github.com/XPIRe-team/XPIRe/archive/master.zip

*On OSX only*, execute the following commands:
 
```
sudo port install gcc48
sudo port select gcc mp-gcc48
sudo port install clang-3.6
sudo port select clang mp-clang-3.6
```

Then execute the following commands to compile everything (boost, gmp, mpfr, create essential files, build client and server):
```
$ cd XPIRe
$ make
```

When this is done you can check that the server/client work correctly with the following commands:
```
$ cd server 
$ ./check-correctness
```

The first test should be pretty long (to build initial caches) and then a set of tests should be display CORRECT or "Skipping test...". If you get INCORRECT tests then something went wrong ...


Usage:
======

XPIRe is composed of a server and a client. Both must be started on their respective directories. Thus to start the server execute:
```
$ cd server
$ ./build/PIRServer
```
And to start the client execute (on a different terminal):
```
$ cd client
$ ./build/PIRClient
```
By default the client tries to reach a local server but a given IP address and port can be specified, use --help to get help on the different options for distant connections.

If run without options the PIR server will look for files in a directory db inside the server directory and consider each file is a database element. The client will present a catalog of the files and ask the user to choose a file. When this is done the client will run an optimizer to decide which are the best cryptographic and PIR parameters to retrieve the file. Then he will send an encrypted PIR Query (i.e. a query that the server will mix with the database without understanding which element it allows to retrieve) to the server. The server then computes an encrypted PIR reply and sends it to the client. Finally, the client will decrypt this reply and store the resulting file in the reception directory inside the client directory.


Available options for the server (PIRServer command):
=====================================================

`-h,--help`                         
Print a help message with the different options.

`-z, --driven`    
Server-driven mode. This mode is to be used when multiple clients will connect to the server with the same cryptographic and PIR parameters. This allows the server to import the database into RAM and to perform precomputations over the database for the first client which **significantly increases the performance for the following clients if LWE-based cryptography is used**. The first client will ask for a given configuration (depending on its optimizer and on the command-line constraints given to the client). After this configuration client, the server will tell the following clients that he is in server-driven mode and that the configuration is imposed. The configuration given by the first client is stored in file arg or in exp/PIRParams.cfg if arg is not specified for further usage (see -L option).

`-L, --load_file arg`     
Load cryptographic and PIR parameters from arg file. Currently unavailable (see issues).

`-s, --split_file arg (=1)`    
Only use first file in db directory and split it in arg database elements. This allows to have a large database with many fixed size elements (e.g. bits, bytes, 24-bit depth points) into a single file which is much more efficient from a file-system point of view than having many small files. Building databases from a single file with more complex approaches (e.g. csv, or sqlite files) would be a great feature to add to XPIRe.

`-p, --port arg (=1234)`    
Port used by the server to listen to incoming connections, by default 1234.

`--db-generator`                        
Generate a fake database with random elements instead of reading it from a directory. This is useful for performance tests. It allows to deal with arbitrary databases without having to build them on the file-system and to evaluate performance costs without considering disk access limitations.

`-n, --db-generator-files arg (=10)`    
Number of files for the virtual database provided by the DB generator.

`-l [ --db-generator-filesize ] arg (=12800000)`    
Filesize in bytes for the files in the virtual database provided by the DB generator.

`--no-pipeline`                         
No pipeline mode. In this mode the server executes each task separately (getting the PIR Query, computing the reply, sending it). Only useful to measure the performance of each step separately.


Available options for the client (PIRClient command):
=====================================================

`-h, --help`                     
Display a help message.

`-i, --serverip arg (=127.0.0.1)`    
Define the IP address at which the client will try to contact the PIRServer.

`-p [ --port ] arg (=1234)`          
Define the port at which the client will try to contact the PIRServer.

`-c, --autochoice`               
Don't display the catalog of database elements and automatically choose the first element without waiting for user input.

`--dry-run`                          
Enable dry-run mode. In this mode the client does not send a PIR Query. It runs the optimizer taking into account the command-line options and outputs the best parameters for each cryptosystem (currently NoCryptography, Paillier and LWE) with details on the costs evaluated for each phase (query generation, query sending, reply generation, reply sending, reply decryption). If a server is available it interacts with it to set the parameters: client-server throughput and server-client throughput. It also requests from the server the performance cache to evaluate how fast the server can process the database for each possible set of cryptographic parameters. If no server is available it uses default performance measures. The other parameters are set for the default example: a thousand mp3 files over ADSL, aggregation disabled and security k=80. Each of these parameters can be overridden on the command line.

`--verbose-optim`                    
Ask the optimizer to be more verbose on the intermediate choices and evaluations (as much output as in the dry-run mode).

`--dont-write`                       
Don't write the result to a file. For testing purposes, it still will process the reply (decryption of the whole answer).

`-f, --file  arg`     
Use a config file to test different optimizations in dry-run mode (see exp/sample.conf). Must be used with the --dry-run option or it is ignored. 


Available options for the optimizer (through PIRClient command):
================================================================

`-n, --file-nbr arg`    
Used in dry-run mode only: Override the default number of database elements.

`-l, --file-size ] arg`    
Used in dry-run mode only: Override the default database element size (in bits).

`-u, --upload arg`     
Force client upload speed in bits/s (bandwidth test will be skipped). This is valid in dry-run or normal mode (e.g. if a user does not want to use more than a given amount of his bandwidth).

`-d, --download arg`       
Force client download speed in bits/s (bandwidth test will be skipped). This is valid in dry-run or normal mode.

`-r, --crypto-params arg`     
Limit with a regular expression arg to a subset of the possible cryptographic parameters. Parameters are dependent on each cryptographic system:
  * NoCryptography if a trivial full database download is to be done after which PIRClient stores only the element the user is interested in.
  * Paillier:A:B:C if Paillier's cryptosystem is to be used with A security bits, a plaintext modulus of B bits and a ciphertext modulus of C bits.
  * LWE:A:B:C if LWE is to be used with A security bits, polynomials of degree B and polynomial coefficients of C bits.
For example it is possible to force just the cryptosystem with NoCryptography.\* or LWE.\*, or ask for a specific parameter set like Paillier:80:1024:2048. Specifying the security with this option is tricky as it must match exactly so better use -k for this purpose.

`-k, --security arg (=80)`     
Minimum security bits required for a set of cryptographic parameters to be considered by the optimizer.

`--dmin arg (=1)`                  
Min dimension value considered by the optimizer. Dimension is also called recursion in the literature. It is done trivially (see the scientific paper) and thus for dimension d query size is proportional to d n^{1/d} and reply size is exponential in d. For databases with many small elements a d>1 can give the best results, but only in exceptional situations having d>4 is interesting.

`--dmax arg (=4)`             
Max dimension value considered by the optimizer.

`-a, --alphaMax  arg (=0)`       
Max aggregation value to test (1 = no aggregation, 0 = no limit). It is sometimes interesting to aggregate a database with many small elements into a database with fewer but larger aggregated elements (e.g. if database elements are one bit long). This value forces the optimizer to respect a maximum value for aggregation, 1 meaning that elements cannot be aggregated.

`-x, --fitness arg (=1)`    
Set fitness method to: 
0=SUM Sum of the times on each task
1=MAX Max of server times + Max of client times
2=CLOUD Dollars in a cloud model (see source code)
This sets the target function of the optimizer. When studying the different parameters the optimizer will choose the one that minimizes this function. 0 corresponds to minimizing the resources spent, 1 to minimizing the round-trip time (given that server operations have are pipelined and client operations are also, independently, pipelined), 2 corresponds to minimizing the cost by associating CPU cycles and bits transmitted to money using a cloud computing model.


Contributors:
=============

This project has been imported to GitHub once mature enough by Carlos Aguilar Melchor which erased all the possible "blames" attributing to each contributor each line of code. In order to give a fair idea of the line contributions per author, the following command:
```
git ls-tree --name-only -z -r HEAD|egrep -z -Z -E '\.(cc|h|cpp|hpp|c|txt)$'  |egrep -z -Z -E -v dependencies|egrep -z -Z -E -v boost|xargs -0 -n1 git blame --line-porcelain|grep "^author "|sort|uniq -c|sort -nr
```
which counts the lines by contributor for the initial commit put into GitHub (removing the dependencies directory which does not correspond to code we developed) gave the following output (aggregating aliases) just before transfering the project to GitHub:
```
7655 author Carlos Aguilar-Melchor
5693 author Joris Barrier
1153 author Marc-Olivier Killijian
```
Besides this line counting, the roles were distributed as follows:    
Carlos Aguilar-Melchor (Associate Professor): co-supervisor, developer    
Joris Barrier (PhD student): developer    
Marc-Olivier Killijian (Researcher): co-supervisor, developer    

Affiliations:
=============
Carlos Aguilar Melchor has been during almost all this project affiliated to the XLIM laboratory in Limoges, he is currently at IRIT laboratory in Toulouse. Joris Barrier and Marc-Olivier Killijian are affiliated to the LAAS-CNRS laboratory in Toulouse. 

The contributors thank Laurent Fousse for his help on evaluating the performance of NTRU.
 
