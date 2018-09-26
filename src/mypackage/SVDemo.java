package com.mypackage;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.ops.CommonOps;

public class SVDemo {
	
	/**
     * <p>
     * SVD decomposition on M and solution of the linear system x=[M I A]
     * </p>
     *
     * @param input1.txt: matrix M 
     * @param input2.txt: vector c
	 * @throws FileNotFoundException, RuntimeException 
	 * @return vector solution of the linear system
	 */
	public static void main(String[] args) throws FileNotFoundException, RuntimeException {
		
		//read matrix M, vector c from i
		double [][]M = readInput(args[0]);
		double [][]C = readInput(args[1]);
		
		//Create DenseMatrix M and c
		DenseMatrix64F m = new DenseMatrix64F(M);
		DenseMatrix64F c = new DenseMatrix64F(C); 
		
		//Print DenseMatrix M and c
	    System.out.println("Printing input matrix M");
	    System.out.println(m.toString());
	    System.out.println("--------------------------------");
	    System.out.println("Printing input vector C");
	    System.out.println(c.toString());
	    System.out.println("--------------------------------");
		
		//SVD
	    SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(m.numRows,m.numCols,true,true,false);	

        if( !svd.decompose(m) )
            throw new RuntimeException("Decosition failed");

        //Get U, S, V
        DenseMatrix64F U = svd.getU(null,false);
        DenseMatrix64F S = svd.getW(null);
        DenseMatrix64F V = svd.getV(null,false);
	    System.out.println("Printing plain U");
	    System.out.println(U.toString());
	    System.out.println("--------------------------------");
	    System.out.println("Printing plain S");
	    System.out.println(S.toString());
	    System.out.println("--------------------------------");
	    System.out.println("Printing plain V");
	    System.out.println(V.toString());
	    System.out.println("--------------------------------");
	    
	    //TEST DECOMPOSITION M = U * S * VT
	    DenseMatrix64F US = new DenseMatrix64F(U.numRows,S.numCols);
	    CommonOps.mult(U, S, US);
	    DenseMatrix64F VT = new DenseMatrix64F(V); 
	   	CommonOps.transpose(VT);
	    DenseMatrix64F USVT = new DenseMatrix64F(US.numRows,VT.numCols);
	    CommonOps.mult(US, VT, USVT);
	    System.out.println("Proving M= U *S * VT");
	    System.out.println(USVT.toString());
	    System.out.println("--------------------------------");
	    
        //Transpose of U
	    DenseMatrix64F UT = CommonOps.transpose(U, null);
	    
	    //Inverse of S, S-1
	    DenseMatrix64F ST = CommonOps.transpose(S, null);
	    int diagLenght = 0;
	    if(S.numRows>S.numCols)
	    	diagLenght=S.numCols;
	    else if(S.numRows<S.numCols)
	    	diagLenght=S.numRows;
	    DenseMatrix64F Sdiag = new DenseMatrix64F(1,diagLenght);
	    CommonOps.extractDiag(ST, Sdiag);
	    double[] diagVal= Sdiag.getData();
	    for(int z=0; z<diagVal.length; z++) {
	    	if(diagVal[z]==0){
	    		diagVal[z]=0;	
	    	}
	    	else{
	    	diagVal[z]=1.0/diagVal[z];
	    	}
	    }
	    Sdiag.setData(diagVal);
	    DenseMatrix64F SI = CommonOps.diagR(ST.numRows, ST.numCols, Sdiag.getData());
	    System.out.println("Printing SI");
	    System.out.println(SI.toString());
	    System.out.println("--------------------------------");
	    
	    //Proving S * SI = I
	    DenseMatrix64F diagonalTest1 = new DenseMatrix64F(S.numRows, SI.numCols);
	    CommonOps.mult(S, SI, diagonalTest1);
	    System.out.println("Proving S * SI = I");
	    System.out.println(diagonalTest1.toString());
	    System.out.println("--------------------------------");
	    
	    //Proving SI * S = I
	    DenseMatrix64F diagonalTest2 = new DenseMatrix64F(SI.numRows, S.numCols);
	    CommonOps.mult(SI, S, diagonalTest2);
	    System.out.println("Proving SI * S = I");
	    System.out.println(diagonalTest2.toString());
	    System.out.println("--------------------------------");
	    
	    // Get M-1 = V * SI * UT
	    DenseMatrix64F VSI = new DenseMatrix64F(V.numRows,SI.numCols);
	    CommonOps.mult(V, SI, VSI);	    	    	
	    DenseMatrix64F VSIUT = new DenseMatrix64F(VSI.numRows,UT.numCols);
	    CommonOps.mult(VSI, UT, VSIUT);
	    System.out.println("Printing M-1 = V * SI * UT");
	    System.out.println(VSIUT.toString());
	    System.out.println("--------------------------------");
	    
	    //  m * M-1 = I
	    DenseMatrix64F matrixTest = new DenseMatrix64F(m.numCols, VSIUT.numRows);
	    DenseMatrix64F mTranspose = new DenseMatrix64F(m); 
	   	CommonOps.transpose(mTranspose);
	    DenseMatrix64F mInvTranspose = new DenseMatrix64F(VSIUT);
	    CommonOps.transpose(mInvTranspose);
	    CommonOps.mult(mTranspose, mInvTranspose, matrixTest);
	    System.out.println("Proving I = m * M-1");
	    System.out.println(matrixTest.toString());
	    System.out.println("--------------------------------");
	    
	    // M-1 * m = I
	    DenseMatrix64F matrixTest2 = new DenseMatrix64F(VSIUT.numRows, m.numCols);
	    CommonOps.mult(VSIUT, m, matrixTest2);
	    System.out.println("Proving I = M-1 * m");
	    System.out.println(matrixTest2.toString());
	    System.out.println("--------------------------------");
	    
	    // Transposing c
	    CommonOps.transpose(c);
	    
	    //x = M-1 * c
	    DenseMatrix64F X = new DenseMatrix64F(VSIUT.numRows, c.numCols);
	    CommonOps.mult(VSIUT, c, X);
	    System.out.println("Printing OUTPUT x of: x = M-1 * c");
	    System.out.println(X.toString());	   
	    System.out.println("--------------------------------");
	}	
	
	private static double[][] readInput (String nameFile) throws FileNotFoundException{
		// read in the data
		Scanner input = new Scanner (new File(nameFile));
		int m=0, n=0;
		// getting matrix dimensions
		if (input.hasNext ()){
			m=input.nextInt();
		}
		if (input.hasNext()){
			n=input.nextInt();
		}
		//creating matrix
		double [][] matrix = new double[m][n];
		for(int x=0; x<m; x++){
			for(int y=0; y<n; y++) {
				if(input.hasNext()){
					matrix[x][y]=input.nextDouble();
				}
			}
		}
		input.close();
		return matrix;
	}
}