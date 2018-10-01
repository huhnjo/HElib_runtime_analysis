#include <iostream>
#include <math.h> 
#include <cmath>
#include<vector>
using namespace std;
#include "FHE.h"
#include "EncryptedArray.h"
#include "EvalMap.h"
#include "powerful.h"
#include <NTL/lzz_pXFactoring.h>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <time.h>
#include "MeasureTime.cpp"

#include "Function_Library_Enc.cpp"

int gcd(int a, int b){
  if (a == 0){
    return b;
  }
  return gcd(b%a, a);
}

// A simple method to evaluate Euler Totient Function
 int phi( int n){
  int res = 1;
  for (int i=2; i < n; ++i)
    if (gcd(i, n) == 1){
      res++;
    }
  return res;
}


string testMult(int R, int x, int y, long m, long phi_m, long d, long mvec1, long mvec2, long gens1, long gens2, long ords1, long ords2, long c_m ){

  // string out;
  std::stringstream out;
  
// ******************** This is setup for the encryption scheme. Encoding stuff starts in line 89. ******************	
// ******************** KeyGen ******************
  //long m=1705;
  //long m = 4369;
  long phim = phi_m;
	long p=2, r=1;
	
	long c=2;
	
	long L=30;
	
	bool cons=1;
	long w=64; 
	//	long d=16;
	ZZX G;
	long dunno=c_m; //The last column in the mValues table in the bootstrapping test - don't know what it does, needed in line 65. 

  Vec<long> mvec;
  vector<long> gens;
  vector<long> ords;
  
  // long phim = phi(m);
  cout << "M=" << m << " Phi(m)=" <<  phim << " d= " << d << endl;
  assert(GCD(p, m) == 1);

  append(mvec, long(mvec1));
  append(mvec, long(mvec2));
  
  gens.push_back(long(gens1));
  gens.push_back(long(gens2));
  
  ords.push_back(long(ords1));
  ords.push_back(long(ords2));

  
  /*append(mvec, long(11));
  append(mvec, long(155));
  
  gens.push_back(long(156));
  gens.push_back(long(936));
  ords.push_back(long(10));
  ords.push_back(long(6));*/

  setTimersOn();
  setDryRun(false);
  double t = -GetTime();
  
  FHEcontext context(m, p, r, gens, ords  );
	context.bitsPerLevel = 23;
	buildModChain(context, L, c);
	context.makeBootstrappable(mvec, /*t=*/0, cons);
	t += GetTime();

	
  long nPrimes = context.numPrimes();
  IndexSet allPrimes(0,nPrimes-1);
  double bitsize = context.logOfProduct(allPrimes)/log(2.0);
  cout << "  "<<nPrimes<<" primes in chain, total bitsize="
       << ceil(bitsize) << ", secparam="
       << (7.2*phim/bitsize -110) << endl;
  //long k = (7.2*phim/bitsize -110);
  //cout << k << endl;
       
  long p2r = context.alMod.getPPowR();
  context.zMStar.set_cM(dunno/100.0);       
	
	FHESecKey secretKey(context);
	FHEPubKey& publicKey = secretKey;
	
	G = makeIrredPoly(p, d); 
    	
    secretKey.GenSecKey(w);
    addSome1DMatrices(secretKey);
    addFrbMatrices(secretKey);
    secretKey.genRecryptData();    
    
  zz_p::init(p2r);
  zz_pX poly_p = random_zz_pX(context.zMStar.getPhiM());
  PowerfulConversion pConv(context.rcData.p2dConv->getIndexTranslation());
  HyperCube<zz_p> powerful(pConv.getShortSig());
  pConv.polyToPowerful(powerful, poly_p);
  ZZX ptxt_poly = conv<ZZX>(poly_p);
  PolyRed(ptxt_poly, p2r, true); // reduce to the symmetric interval    
    
    
    // ******************** Make plaintexts *********************
    
  int a = x;  // Value one as integer
  int b = y;	// Value two as integer
	
	vector<bool> a1 = IntToTwoCompMin(a);  	// Transform a into its Two's Complement representation
	vector<bool> b1 = IntToTwoCompMin(b);		// Transform b into its Two's Complement representation 
	
	int na=a1.size();
	int nb=b1.size();
	
	int a2 = TwoCompToInt(a1);			// Convert values back to check that conversion worked correctly
	int b2 = TwoCompToInt(b1);
	
	cout<<a<<", "<<a2<<'\n';
	cout<<b<<", "<<b2<<'\n';
    
    // ******************** Encryption ******************
	EncryptedArray ea(context, G);    
    
    // This takes the boolean array and converts in to the right datatype for plaintexts.
    
	vector<long> m1[na];
	for(int i=0;i<na;i++){
		m1[i].push_back(a1[i]);
		while(int(m1[i].size())<ea.size()){
			m1[i].push_back(0);
		}
	}
	
	vector<long> m2[nb];
	for(int i=0;i<nb;i++){
		m2[i].push_back(b1[i]);
		while(int(m2[i].size())<ea.size()){
			m2[i].push_back(0);
		}
	}
	
	// The following two are just the bits 0 and 1 encrypted, respectively (we sometimes need them)
	vector<long> mbit0;
	mbit0.push_back(0);
	while(int(mbit0.size())<ea.size()){
		mbit0.push_back(0);
	}	
	
	vector<long> mbit1;
	mbit1.push_back(1);
	while(int(mbit1.size())<ea.size()){
		mbit1.push_back(0);
	}		
    	
    cout<<"Message Generation done \n";
    
    // Encrypt the plaintexts:
    vector<Ctxt> ct1;
    for(int i=0;i<na;i++){
    	Ctxt foo(publicKey);
    	ea.encrypt(foo, publicKey, m1[i]);
    	ct1.push_back(foo);
	}


    vector<Ctxt> ct2;
    for(int i=0;i<nb;i++){
    	Ctxt foo(publicKey);
    	ea.encrypt(foo, publicKey, m2[i]);
    	ct2.push_back(foo);
	}     
    
    Ctxt cbit0(publicKey);
	ea.encrypt(cbit0, publicKey, mbit0);

    Ctxt cbit1(publicKey);
	ea.encrypt(cbit1, publicKey, mbit1);	
    	
    cout<<"Enc done \n"; 
    
    
    
    
    // ******************** Computation ******************
    MeasureTime timer;
    
    //vector<Ctxt> ctsum = TCEncAdd(ct1,ct2); // Add ct1 and ct2, no bootstraaping.
   // With bootstrapping: ctsum = TCEncBootAdd(publicKey, ct1,ct2)    
    vector<Ctxt> ctprod = ct1;

    for(int i = 0; i < R; i++){
      timer.startMeasuring();
    
      ctprod = FancyMult(publicKey, ctprod, ct2, cbit0, cbit1); // Multiply ct1 and ct2, always with bootstrapping.

      long elapsed_time = timer.endMeasuring();
      
      out  << elapsed_time << endl;
      
      cout << "TIME COMP ITERATION " << i <<  "  -  " << elapsed_time << endl;

    }
    // ctprod = FancyMult(publicKey, ct1ctprod, ct2, cbit0, cbit1); // Multiply ct1 and ct2, always with bootstrapping.
   

     
    cout<<"Computation done \n";
        
    
    // ******************** Decrypt ******************

 	vector<long> foo;

	vector<bool> resadd;
	vector<bool> resmult;
	
	
	// Decrypt result of addition
	/*    for(int i=0;i<int(ctsum.size());i++){
    	ea.decrypt(ctsum[i], secretKey, foo);
    	resadd.push_back(bool(foo[0]));
	}*/

	// Decrypt result of multiplication
    for(int i=0;i<int(ctprod.size());i++){
    	ea.decrypt(ctprod[i], secretKey, foo);
    	resmult.push_back(bool(foo[0]));
	}

	for(int i=0;i<resmult.size();i++){
		cout<<resmult[i];
	}
	cout<<endl;

	// Decode back from Boolean vector to int
	//int decadd = TwoCompToInt(resadd); 
	int decmult = TwoCompToInt(resmult);
	
	//int rightadd = a+b;
	int rightmult = a*b;
	
	// Check if results are correct
	//	cout<<"a+b: "<<decadd<<" = "<<rightadd<<endl;
	cout<<"a*b: "<<decmult<<" = "<<rightmult<<endl;

  
	return out.str();
}



int main(){

  long p = 2;
  

   long table[][14] = { 
//{  p, phi(m),    m, d,  m1,  m2,  m3,   g1,    g2,   g3, ord1,ord2,ord3,c_m}
  {  2,   600,  1023, 10, 11,  93,  0,   838,   584,    0, 10,  6,   0, 100}, 
  {  2,  1200,  1705, 20, 11, 155,  0,   156,   936,    0, 10,  6,   0, 100}, 
  {  2,  4096,  4369, 16, 17, 257,  0,   258,  4115,    0, 16,-16,   0, 100}, 
  {  2, 12800, 17425, 40, 41, 425,  0,  5951,  8078,    0, 40, -8,   0, 100}, 
  //{  2, 15004, 15709, 22, 23, 683,  0,  4099, 13663,    0, 22, 31,   0, 100}, 
  //{  2, 18000, 18631, 25, 31, 601,  0, 15627,  1334,    0, 30, 24,   0,  50}, 
  {  2, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0, 100}, 
  // {  2, 21168, 27305, 28, 43, 635,  0, 10796, 26059,    0, 42, 18,   0, 100}, 
  //{  2, 24000, 31775, 20, 41, 775,  0,  6976, 24806,    0, 40, 30,   0, 100}, 
  //{  2, 26400, 27311, 55, 31, 881,  0, 21145,  1830,    0, 30, 16,   0, 100}, 
  {  2, 31104, 35113, 36, 37, 949,  0, 16134,  8548,    0, 36, 24,   0, 400}, 
  {  2, 34848, 45655, 44, 23, 1985, 0, 33746, 27831,    0, 22, 36,   0, 100}, 
  //{  2, 42336, 42799, 21, 127, 337, 0, 25276, 40133,    0,126, 16,   0,  20}, 
  //{  2, 45360, 46063, 45, 73, 631,  0, 35337, 20222,    0, 72, 14,   0, 100}, 
  //{  2, 49500, 49981, 30, 151, 331, 0,  6952, 28540,    0,150, 11,   0, 100}, 
  {  2, 54000, 55831, 25, 31, 1801, 0, 19812, 50593,    0, 30, 72,   0, 100}, 
  //{  2, 60016, 60787, 22, 89, 683,  0,  2050, 58741,    0, 88, 31,   0, 200}
   };



   int j = 0;
   ofstream runtimeResults;
   runtimeResults.open ("polynomialConstantRes.txt");

   int a = 4;
   int b = 2;
   int R = 10;

     string out = "";
     /*
   for(int i = 0; i < R; i++){
     b = pow(2, i+1);
     out = testMult(1, a, b, table[j][2], table[j][1], table[j][3], table[j][4], table[j][5], table[j][7], table[j][8], table[j][10], table[j][11], table[j][13] );
     runtimeResults << a << "," << b << "," << out << endl;
   }
   
   a=8;   
   for(int i = 0; i < R; i++){
     b = pow(2, i+1);
     out = testMult(1, a, b, table[j][2], table[j][1], table[j][3], table[j][4], table[j][5], table[j][7], table[j][8], table[j][10], table[j][11], table[j][13] );
     runtimeResults << a << "," << b << "," << out << endl;
     }*/

   a=16;
   for(int i = 0; i < R; i++){
     b = pow(2, i+1);
     out = testMult(1, a, b, table[j][2], table[j][1], table[j][3], table[j][4], table[j][5], table[j][7], table[j][8], table[j][10], table[j][11], table[j][13] );
     runtimeResults << a << "," << b << "," << out << endl;
   }
   /*
   a=128;
   for(int i = 0; i < R; i++){
     b = pow(2, i+1);
     out = testMult(1, a, b, table[j][2], table[j][1], table[j][3], table[j][4], table[j][5], table[j][7], table[j][8], table[j][10], table[j][11], table[j][13] );
     runtimeResults << a << "," << b << "," << out << endl;
     }*/
   
   runtimeResults.close();
   
   return 0;
}
