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


long testMult(int x, int y, long R ){

  
// ******************** This is setup for the encryption scheme. Encoding stuff starts in line 89. ******************	
// ******************** KeyGen ******************
	long m=1705, p=2, r=1;
	long c=3;
	long L=40;
	long B=23;
//	long N=0;
//	long skHwt=5;
	bool cons=1;
//	long nthreads=4;
//	long seed=0;	
	long w=64; 
	long d=20;
	ZZX G;
	long dunno=100; //The last column in the mValues table in the bootstrapping test - don't know what it does, needed in line 65. 

  Vec<long> mvec;
  vector<long> gens;
  vector<long> ords;
  
  long phim = phi(m);
  cout << "M=" << m << " Phi(m)=" <<  phim << endl;
  assert(GCD(p, m) == 1);

  append(mvec, long(11));
  append(mvec, long(155));
  gens.push_back(long(156));
  gens.push_back(long(936));
  ords.push_back(long(10));
  ords.push_back(long(6));

  setTimersOn();
  setDryRun(false);
  double t = -GetTime();  
	
	FHEcontext context(m, p, r, gens, ords);
	context.bitsPerLevel = B;
	buildModChain(context, L, c);
	context.makeBootstrappable(mvec, /*t=*/0, cons);
	t += GetTime();
	
  long nPrimes = context.numPrimes();
  IndexSet allPrimes(0,nPrimes-1);
  double bitsize = context.logOfProduct(allPrimes)/log(2.0);
  cout << "  "<<nPrimes<<" primes in chain, total bitsize="
       << ceil(bitsize) << ", secparam="
       << (7.2*phim/bitsize -110) << endl;	
       
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
    timer.startMeasuring();
    vector<Ctxt> ctprod;
    for(int i = 0; i < R; i++){
      ctprod = FancyMult(publicKey, ct1, ct2, cbit0, cbit1); // Multiply ct1 and ct2, always with bootstrapping.
    }
    long elapsed_time = timer.endMeasuring();
     
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
	
	int rightadd = a+b;
	int rightmult = a*b;
	
	// Check if results are correct
	//	cout<<"a+b: "<<decadd<<" = "<<rightadd<<endl;
	cout<<"a*b: "<<decmult<<" = "<<rightmult<<endl;

  
	return elapsed_time;
}


long testPoly(int x, int y ){

  
// ******************** This is setup for the encryption scheme. Encoding stuff starts in line 89. ******************	
// ******************** KeyGen ******************
	long m=1705, p=2, r=1;
	long c=3;
	long L=40;
	long B=23;
//	long N=0;
//	long skHwt=5;
	bool cons=1;
//	long nthreads=4;
//	long seed=0;	
	long w=64; 
	long d=20;
	ZZX G;
	long dunno=100; //The last column in the mValues table in the bootstrapping test - don't know what it does, needed in line 65. 

  Vec<long> mvec;
  vector<long> gens;
  vector<long> ords;
  
  long phim = phi(m);
  cout << "M=" << m << " Phi(m)=" <<  phim << endl;
  assert(GCD(p, m) == 1);

  append(mvec, long(11));
  append(mvec, long(155));
  gens.push_back(long(156));
  gens.push_back(long(936));
  ords.push_back(long(10));
  ords.push_back(long(6));

  setTimersOn();
  setDryRun(false);
  double t = -GetTime();  
	
	FHEcontext context(m, p, r, gens, ords);
	context.bitsPerLevel = B;
	buildModChain(context, L, c);
	context.makeBootstrappable(mvec, /*t=*/0, cons);
	t += GetTime();
	
  long nPrimes = context.numPrimes();
  IndexSet allPrimes(0,nPrimes-1);
  double bitsize = context.logOfProduct(allPrimes)/log(2.0);
  cout << "  "<<nPrimes<<" primes in chain, total bitsize="
       << ceil(bitsize) << ", secparam="
       << (7.2*phim/bitsize -110) << endl;	
       
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
  //  int b = y;	// Value two as integer
	
	vector<bool> a1 = IntToTwoCompMin(a);  	// Transform a into its Two's Complement representation
	//	vector<bool> b1 = IntToTwoCompMin(b);		// Transform b into its Two's Complement representation 
	
	int na=a1.size();
	//	int nb=b1.size();
	
	int a2 = TwoCompToInt(a1);			// Convert values back to check that conversion worked correctly
	//	int b2 = TwoCompToInt(b1);
	
	cout<<a<<", "<<a2<<'\n';
	//	cout<<b<<", "<<b2<<'\n';
    
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
	
	/*	vector<long> m2[nb];
	for(int i=0;i<nb;i++){
		m2[i].push_back(b1[i]);
		while(int(m2[i].size())<ea.size()){
			m2[i].push_back(0);
		}
		}*/
	
	// The following two are just the bits 0 and 1 encrypted, respectively (we sometimes need them)
	vector<long> mbit0;
	mbit0.push_back(0);
	while(int(mbit0.size())<ea.size()){
		mbit0.push_back(0);
	}	
	
	/*	vector<long> mbit1;
	mbit1.push_back(1);
	while(int(mbit1.size())<ea.size()){
		mbit1.push_back(0);
		}*/		
    	
    cout<<"Message Generation done \n";
    
    // Encrypt the plaintexts:
    vector<Ctxt> ct1;
    for(int i=0;i<na;i++){
    	Ctxt foo(publicKey);
    	ea.encrypt(foo, publicKey, m1[i]);
    	ct1.push_back(foo);
	}


    /*    vector<Ctxt> ct2;
    for(int i=0;i<nb;i++){
    	Ctxt foo(publicKey);
    	ea.encrypt(foo, publicKey, m2[i]);
    	ct2.push_back(foo);
	}     */
    
    Ctxt cbit0(publicKey);
	ea.encrypt(cbit0, publicKey, mbit0);

	/*    Ctxt cbit1(publicKey);
	      ea.encrypt(cbit1, publicKey, mbit1);	*/
    	
    cout<<"Enc done \n"; 
    
    
    
    
    // ******************** Computation ******************
    MeasureTime timer;
    
    //vector<Ctxt> ctsum = TCEncAdd(ct1,ct2); // Add ct1 and ct2, no bootstraaping.
   // With bootstrapping: ctsum = TCEncBootAdd(publicKey, ct1,ct2)
    timer.startMeasuring();
    vector<Ctxt> ctprod;

    ctprod = FancyMult(publicKey, ct1, ct1, cbit0, cbit0); // Multiply ct1 and ct2, always with bootstrapping.

    long elapsed_time = timer.endMeasuring();
     
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
	
	//	int rightadd = a+b;
	//int rightmult = a*;
	
	// Check if results are correct
	//	cout<<"a+b: "<<decadd<<" = "<<rightadd<<endl;
	cout<<"Result = "<<decmult << endl;

  
	return elapsed_time;
}


int main(){

  /*
// multiply by constant.
  ofstream runtimeResults;
  runtimeResults.open ("runtimeResultsConstantMult.txt");

  int a = 2;
  int power_of_a = 0;
  int b = 1;

  do{

    int power = pow(2, power_of_a);

    for(int c = 1; c < 4; c++){
      long time = testMult(power, b, 1);
      runtimeResults << "a=" << power <<  "," << "b=" << b << "," << time << endl;     
    }
    
    power_of_a++;
  }while(power_of_a <= 30);

  runtimeResults.close();
  
return 0;
  */

  /******* check multiplication regression behaviour*****
  ofstream runtimeResults;
  runtimeResults.open ("runtimeResultsMultRegression.txt");

  int fixedPoints [10] = {pow(2,1), pow(2,3), pow(2,6), pow(2,11), pow(2,15), pow(2,18),  pow(2,22), pow(2,25), pow(2,27), pow(2,30)};
  int a = 128;

  for(int i = 0; i < 10; i++){
    long time = testMult(a, fixedPoints[i], 1 );
    runtimeResults << "a=" << a << "," << "b=" << fixedPoints[i] << "," << time << endl;
  }

  runtimeResults.close();

  return 0;
  */


  /** test poly*/
   ofstream runtimeResults;
  runtimeResults.open ("polyResult.txt");

  
  int a = 4;
  int exponent = 10;
  
  
  for(int i = 0; i < 10; i++){

    int b = 4;

    long time = testMult(a,b, 1);

    runtimeResults << "a=" << a << "," << "b=" << b << "," << time << endl;

    int result = a*b;
    
    
    a = result;

  }
  
}
