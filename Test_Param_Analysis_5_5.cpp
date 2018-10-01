#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include <NTL/lzz_pXFactoring.h>
#include "MeasureTime.cpp"
#include "Ctxt.cpp"
#include <iostream>
#include <fstream>
#include <cassert>

class Test_Multiplication {

public:
  Test_Multiplication(){

  }

  MeasureTime multiplicationTime;
  long elapsed_time;
  
	

  void  TestCase(long R, long p, long r, long d, long c, long k, long w, long L, long m)
	{

	  cout << "--------------------------" << endl;

	  // cout << "FHE SIZE:= " << FHE_p2Size << endl;
	  
	  std::cout     << "R="  << R    //R = Number of Rounds
			<< ", p=" << p    //p = Plaintextbase
			<< ", r=" << r    //r = lifting
			<< ", d=" << d    //d = 
			<< ", c=" << c
			<< ", k=" << k
			<< ", w=" << w    //Hamming - weight
			<< ", L=" << L
			<< ", m=" << m
			<< endl;



	   vector<long> gens1, ords1;;	  
	   FHEcontext context(m, p, r, gens1, ords1);
	   buildModChain(context, L, c);

	   // ************ calculate G ***********	   
	    ZZX G;
	    if (d == 0)
	      G = context.alMod.getFactorsOverZZ()[0];
	    else
	      G = makeIrredPoly(p, d); 
	    cout << "G = " << G << endl;
	    cout << "Security Level = " << ceil(context.securityLevel()) << endl;
	    // ************ end G ****************

	    FHESecKey secretKey(context);
	    const FHEPubKey& publicKey = secretKey;
	    secretKey.GenSecKey(w);                             //w=Hamming weight
	    addSome1DMatrices(secretKey);

	    EncryptedArray ea(context, G);
	    long nslots = ea.size();
	    cout << "Slots=" << nslots << endl;

	    NewPlaintextArray p0(ea);
	    NewPlaintextArray p1(ea);

	    random(ea, p0);
	    random(ea, p1);

	    cout << p0 << endl;

	    Ctxt c0(publicKey), c1(publicKey);

	    ea.encrypt(c0, publicKey, p0);
	    ea.encrypt(c1, publicKey, p1);

	   
	    //cout << c1.findBaseLevel() << endl;

	     NewPlaintextArray const1(ea);
	     random(ea, const1);
	     
	     ZZX const1_poly;
	     ea.encode(const1_poly, const1);

	     for(long i = 0; i < R; i++) {
	       mul(ea, p1, p0);
	     }

	     multiplicationTime.startMeasuring();
	     for (long i = 0; i < R; i++) {	     
	       c1.multiplyBy(c0);
	     }
	     elapsed_time = multiplicationTime.endMeasuring();
	     cout << "Elapsed time: " << elapsed_time << " milliseconds" << endl;

	     c0.cleanUp();
	     c1.cleanUp();

	     NewPlaintextArray pp0(ea);
	     NewPlaintextArray pp1(ea);

	     ea.decrypt(c0, secretKey, pp0);
	     ea.decrypt(c1, secretKey, pp1);

	     
	     cout << "Results: " << std::boolalpha << endl;
	     cout << "pp0 equals p0 ? " <<  equals(ea, pp0, p0) <<   endl;
	     cout << "pp1 equals p1 ? " << equals(ea, pp1, p1) << endl;
	     cout << "---------------------------" << endl;
	    
	}

  long calculateSlots( long p, long r, long d, long c, long k, long w, long L, long m ){
           vector<long> gens1, ords1;
	   FHEcontext context(m, p, r, gens1, ords1);
	   buildModChain(context, L, c);

	   // ************ calculate G ***********	   
	    ZZX G;
	    if (d == 0)
	      G = context.alMod.getFactorsOverZZ()[0];
	    else
	      G = makeIrredPoly(p, d); 
	    //  cout << "G = " << G << endl;
	    //	    cout << "Security Level = " << ceil(context.securityLevel()) << endl;
	    // ************ end G ****************

	    FHESecKey secretKey(context);
	    const FHEPubKey& publicKey = secretKey;
	    secretKey.GenSecKey(w);                             //w=Hamming weight
	    addSome1DMatrices(secretKey);

	    EncryptedArray ea(context, G);
	    long nslots = ea.size();
	    return nslots;
  }

  
};

int main()
{
        Test_Multiplication tm;

	/*default values */
	/*	long p = 2;
	long r = 1;
	long d = 1;
	long c = 2;
	long k = 80;
	long w = 64;
	long s = 0;

	
	long R = 10;
	long L = 0;
	
	if (L==0) { 
	    L = R+2;
	}
	
	
	long chosen_m = 0;
	long m = FindM(k, L, c, p, d, s, chosen_m, false);
	
	
	for(long i = R; i > 0; i--){
	  L = i + 2;
	  tm.TestCase(i /*=R ,p,r,d,c,k,w,L,m);
	}*/


	long p = 2;
	long r = 1;
	long d = 1;
	long c = 2;
	long k = 80;
	long w = 64;
	long s = 0;

	
	long R = 1;
	long L = 0;

	L = R+2;

        long chosen_m = 0; 
	long two_i = 2;
	long noOfSlots = 0;

	ofstream lookupFile;
	lookupFile.open ("LookUpTableForSlots.txt");
	
	
	for(int j = 1; j < 40; j++){
	  L = j+2;
	  
	  for(int i = 0; i < 11; i++){
	    two_i = pow(2,i); // minimum no. Of Slots
	    
	  
	    long m = FindM(k, L, c, p, d, two_i, chosen_m, false);
	    noOfSlots = tm.calculateSlots( p,  r,  d,  c,  k,  w,  L, m);
	    
	     cout << "L=" << L << " | Min. No. of Slots= " << two_i << " | Real No. of Slots=" << noOfSlots << " | Modulus= " << m << endl;
	     lookupFile << L << "," << two_i << "," << noOfSlots << "," << m << endl;
	  }
	  
	}

	lookupFile.close();
	
	
	return 0;
}


