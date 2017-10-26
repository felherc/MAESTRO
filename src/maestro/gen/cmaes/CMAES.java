/**
 * Copyright 2017 Felipe Hernández
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software distributed under the 
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 * express or implied. See the License for the specific language governing permissions and
 * limitations under the License.
 */

package maestro.gen.cmaes;

import java.util.ArrayList;

import maestro.ContVar;
import maestro.DiscVar;
import maestro.gen.Generator;
import maestro.solution.Solution;
import maestro.solution.SolutionRoot;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import probDist.Uniform;
import probDist.multiVar.MultiVarNormal;
import utilities.MatUtil;
import utilities.stat.ContSeries;

/**
 * This generator allows to create new solutions to a problem with discrete and continuous 
 * variables using a sorted base population of solutions and an Covariance Matrix Adaptation
 * Evolution Strategy (CMA-ES) algorithm
 * @author Felipe Hernández
 */
public class CMAES implements Generator
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The identifier of the generator method
	 */
	public final static String ID = "Covariance Matrix Adaptation Evolution Strategy";
	
	/**
	 * The short version of the identifier of the generator method
	 */
	public final static String SHORT_ID = "CMA-ES";
	
	/**
	 * The <code>q</code> attribute printing name
	 */
	public final static String PARAM_Q = "q = ";
	
	/**
	 * The <code>trunc</code> attribute printing name
	 */
	public final static String PARAM_TRUNC = "Truncated percentage = ";
	
	/**
	 * Default value for the <code>q</code> attribute
	 */
	public final static double D_Q = 0.2;
	// Important - sweet spot between 0.1 (10) and 0.2 (30) - higher for more dimensions
	
	/**
	 * Default value for the <code>trunc</code> attribute
	 */
	public final static double D_TRUNC = 0.45;
	// Sweet spot between 0.4 and 0.8; a lot of variation for higher dimensions
	// Interaction with q
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The list of the discrete variables of the problem
	 */
	private ArrayList<DiscVar> discVars;
	
	/**
	 * The list of the continuous variables of the problem
	 */
	private ArrayList<ContVar> contVars;
	
	/**
	 * This variable defines how the solutions in the population are assigned weights. The 
	 * weights are used to compute the multivariate normal distribution to generate the new 
	 * solutions from. <i>q</i> has to be positive and bigger than 0. A value close to <i>0</i> 
	 * means only the most fitted solutions may contribute to the distribution. A large value of 
	 * <i>q</i> means that every solution has the same contribution with little regard for its 
	 * ranking. Use <code>Double.NaN</code> if an exact uniform distribution should be used.
	 */
	private double q;
	
	/**
	 * The current iteration number, starting from 0
	 */
	private int iteration;
	
	/**
	 * The size of the population in the last call. If the current population size is different, a
	 * new process starts by initializing the control variables.
	 */
	private int prevPop;
	
	/**
	 * The percentage of solutions that are to be truncated from selection. The 
	 * <i>(1 - <code>trunc</code>)*n</i> best solutions are selected.
	 */
	private double trunc;
	
	/**
	 * The parameter for weighting between rank-one and rank-<i>miu</i> update
	 */
	private double miuCov;
	
	/**
	 * The number of parents selected to create new solutions
	 */
	private int miu;
	
	/**
	 * The variance effective selection mass
	 */
	private double miuEff;
	
	/**
	 * Learning rate for cumulation for the rank-one update of the covariance matrix. Must be
	 * greater than one.
	 */
	private double cC;
	
	/**
	 * Learning rate for the cumulation for the step size control. Must be greater than one.
	 */
	private double cSigma;
	
	/**
	 * Learning rate for the covariance matrix update. Must be greater than one.
	 */
	private double cCov;
	
	/**
	 * Damping parameter for step size update. Should be approximately one.
	 */
	private double dSigma;
	
	/**
	 * The weights of the solutions in the current population
	 */
	private double[] weights;
	
	/**
	 * Mean value of the search distribution from the last generation
	 */
	private double[] mean;
	
	/**
	 * Evolution path, a sequence of successive (normalized) steps, for the step size control from 
	 * the last generation
	 */
	private double[] pSigma;
	
	/**
	 * Evolution path, a sequence of successive (normalized) steps, for the rank-one update of the
	 * covariance matrix from the last generation
	 */
	private double[] pC;
	
	/**
	 * Step size control from the last generation
	 */
	private double sigma;
	
	/**
	 * Covariance matrix from the last generation
	 */
	private double[][] C;
	
	/**
	 * The inverse of the square root of the covariance matrix (<i>C<sup>-1/2</sup></i>) from the 
	 * last generation
	 */
	private RealMatrix CSqrtInv;
	
	/**
	 * The number of generations that the <code>CSqrtInv</code> has lived without being updated
	 */
	private int cSqrtInvLife;
	
	// --------------------------------------------------------------------------------------------
	// Constructors
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new instance of the Covariance matrix adaptation evolution strategy generator with 
	 * the default parameters
	 */
	public CMAES()
	{
		discVars = new ArrayList<DiscVar>();
		contVars = new ArrayList<ContVar>();
		q = D_Q;
		trunc = D_TRUNC;
	}
	
	/**
	 * Creates a new instance of the Covariance matrix adaptation evolution strategy generator
	 * @param q This variable defines how the solutions in the population are assigned weights. The 
	 * weights are used to compute the multivariate normal distribution to generate the new 
	 * solutions from. <i>q</i> has to be positive and bigger than 0. A value close to <i>0</i> 
	 * means only the most fitted solutions may contribute to the distribution. A large value of 
	 * <i>q</i> means that every solution has the same contribution with little regard for its 
	 * ranking. Use <code>Double.NaN</code> if an exact uniform distribution should be used.
	 * @param trunc The percentage of solutions that are to be truncated from selection. The 
	 * <i>(1 - <code>trunc</code>)*n</i> best solutions are selected.
	 */
	public CMAES(double q, double trunc)
	{
		discVars = new ArrayList<DiscVar>();
		contVars = new ArrayList<ContVar>();
		if(Double.isNaN(q))
			this.q = q;
		else
			this.q = q <= 0 ? 1E-12 : q;
		this.trunc = trunc < 0 ? 0.0 : trunc;
		this.trunc = trunc > 1 ? 1.0 : this.trunc;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	@Override
	public String getId() 
	{
		return ID;
	}

	@Override
	public String getShortId() 
	{
		return SHORT_ID;
	}
	
	/**
	 * @return This variable defines how the solutions in the population are assigned weights. The 
	 * weights are used to compute the multivariate normal distribution to generate the new 
	 * solutions from. <i>q</i> has to be positive and bigger than 0. A value close to <i>0</i> 
	 * means only the most fitted solutions may contribute to the distribution. A large value of 
	 * <i>q</i> means that every solution has the same contribution with little regard for its 
	 * ranking. Use <code>Double.NaN</code> if an exact uniform distribution should be used.
	 */
	public double getQ() 
	{
		return q;
	}

	/**
	 * @param q This variable defines how the solutions in the population are assigned weights. The 
	 * weights are used to compute the multivariate normal distribution to generate the new 
	 * solutions from. <i>q</i> has to be positive and bigger than 0. A value close to <i>0</i> 
	 * means only the most fitted solutions may contribute to the distribution. A large value of 
	 * <i>q</i> means that every solution has the same contribution with little regard for its 
	 * ranking. Use <code>Double.NaN</code> if an exact uniform distribution should be used.
	 */
	public void setQ(double q) 
	{
		if(Double.isNaN(q))
			this.q = q;
		else
			this.q = q <= 0 ? 1E-12 : q;
	}

	/**
	 * @return The percentage of solutions that are to be truncated from selection. The 
	 * <i>(1 - <code>trunc</code>)*n</i> best solutions are selected.
	 */
	public double getTrunc() 
	{
		return trunc;
	}

	/**
	 * @param trunc The percentage of solutions that are to be truncated from selection. The 
	 * <i>(1 - <code>trunc</code>)*n</i> best solutions are selected.
	 */
	public void setTrunc(double trunc) 
	{
		this.trunc = trunc < 0 ? 0.0 : trunc;
		this.trunc = trunc > 1 ? 1.0 : this.trunc;
	}

	@Override
	public String getParamSummary() 
	{
		return PARAM_Q + q + "; " + PARAM_TRUNC + trunc;
	}

	@Override
	public void addDiscVariable(DiscVar variable) 
	{
		discVars.add(variable);
	}

	@Override
	public void addContVariable(ContVar variable) 
	{
		contVars.add(variable);
	}

	@Override
	public void clearVariables() 
	{
		discVars = new ArrayList<DiscVar>();
		contVars = new ArrayList<ContVar>();
	}

	@Override
	public ArrayList<SolutionRoot> generateSolutions(ArrayList<Solution> population, int number)
	{
		int popSize = population.size();
		int discVarCount = discVars.size();
		int contVarCount = contVars.size();
		int n = discVarCount + contVarCount;
		
		if(popSize == 0 || n == 0)
			return null;
				
		// Initialize attributes if necessary
		if(popSize != prevPop)
			initAttributes(n, popSize);
		prevPop = popSize;
		
		// Compute mean
		ArrayList<ContSeries> varMeans = new ArrayList<ContSeries>();
		for(int j = 0 ; j < n ; j++)
			varMeans.add(new ContSeries(true));
		for(int i = 0 ; i < miu ; i++)
		{
			Solution solution = population.get(i);
			for(int j = 0 ; j < n ; j++)
			{
				double value = j < discVarCount ? solution.getDiscValues().get(j)
												: solution.getContValues().get(j - discVarCount);
				varMeans.get(j).addValue(value);
			}
		}
		double[] newMean = new double[n];
		for(int j = 0 ; j < n ; j++)
			newMean[j] = varMeans.get(j).getMean();
		
		// Initialize temporal variables
		boolean validC = true;
		double[] newPSigma = null;
		double[] newPC = null;
		double newSigma = Double.NaN;
		double[][] newC = null;
		double[][] genC = null;
		
		try
		{		
			// Step size control - pSigma
			if(cSqrtInvLife >= Math.max(1, (int)(1/(10.0*(double)n*cCov))))
			{
				RealMatrix Cmat = new Array2DRowRealMatrix(C);
				EigenDecomposition ed = new EigenDecomposition(Cmat);
				RealMatrix Bmat = ed.getV();
				RealMatrix BTmat = ed.getVT();
				double[][] D = ed.getD().getData();
				for(int i = 0 ; i < D.length ; i++)
					D[i][i] = 1.0/Math.sqrt(D[i][i]); // Inverse of the square root of D
				RealMatrix Dmat = new Array2DRowRealMatrix(D);
				CSqrtInv = Bmat.multiply(Dmat.multiply(BTmat));
				cSqrtInvLife = 0;
			}
			cSqrtInvLife++;
			double esc = Math.sqrt(cSigma*(2 - cSigma)*miuEff);
			double[] vM = new double[n];
			for(int i = 0 ; i < n ; i++)
				vM[i] = (newMean[i] - mean[i])/sigma;
			double[] v1 = CSqrtInv.scalarMultiply(esc).operate(vM);
			newPSigma = new double[n];
			for(int i = 0 ; i < n ; i++)
				newPSigma[i] = (1 - cSigma)*pSigma[i] + v1[i];
			
			// Step size control - sigma
			double newPSigmaNorm = new Array2DRowRealMatrix(newPSigma).getNorm();
			double eN0I = Math.sqrt(n)*(1 - 1/(4*n) + 1/(21*n*n));
			newSigma = sigma*Math.exp(cSigma/dSigma*(newPSigmaNorm/eN0I - 1));
			
			// Covariance matrix adaptation - pC
			newPC = new double[n];
			double hSigma = newPSigmaNorm/Math.sqrt(1 - (Math.pow(1 - cSigma, 2*(iteration + 1))))
							< (1.5 + 1/(n - 0.5))*eN0I ? 1 : 0;
			double dHSigma = (1 - hSigma)*cC*(2 - cC);
			for(int i = 0 ; i < n ; i++)
				newPC[i] = (1 - cC)*pC[i] + hSigma*esc*vM[i];
			
			// Covariance matrix adaptation - newC
			newC = new double[n][n];
			double esc1 = 1 - cCov;
			double esc2 = cCov/miuCov;
			double esc3 = cCov*(1 - 1/miuCov);
			double[][] oPPC = MatUtil.outerProduct(newPC, newPC);
			double[][] sumOPDev = new double[n][n];
			for(int parent = 0 ; parent < miu ; parent++)
			{
				double[] dev = new double[n];
				ArrayList<Integer> discValues = population.get(parent).getDiscValues();
				int var = 0;
				for(int i = 0 ; i < discVars.size() ; i++)
				{
					double value = discValues.get(i);
					dev[var] = (value - mean[var])/sigma;
					var++;
				}
				ArrayList<Double> contValues = population.get(parent).getContValues();
				for(int i = 0 ; i < contVars.size() ; i++)
				{
					double value = contValues.get(i);
					dev[var] = (value - mean[var])/sigma;
					var++;
				}
				double[][] oPDev = MatUtil.outerProduct(dev, dev);
				for(int i = 0 ; i < n ; i++)
					for(int j = 0 ; j < n ; j++)
						sumOPDev[i][j] = sumOPDev[i][j] + weights[parent]*oPDev[i][j];
			}
			genC = new double[n][n];
			for(int i = 0 ; i < n ; i++)
			{
				for(int j = 0 ; j < n ; j++)
				{
					newC[i][j] = esc1*C[i][j] + esc2*(oPPC[i][j] + dHSigma*C[i][j]) 
									+ esc3*sumOPDev[i][j];
					genC[i][j] = newSigma*newSigma*newC[i][j];
					if(Double.isNaN(genC[i][j]) || Double.isInfinite(genC[i][j]))
					{
						validC = false;
						break;
					}
				}
				if(!validC)
					break;
			}
		} catch(Exception e) // The new matrices could not be computed
		{
			validC = false;
		}
		
		// Update variables
		if(validC)
		{
			iteration++;
			mean = newMean;
			pSigma = newPSigma;
			pC = newPC;
			sigma = newSigma;
			C = newC;
		}
		else // Re-initialize the matrices
		{
			initAttributes(n, popSize);
			mean = newMean;
			genC = new double[n][n];
			for(int i = 0 ; i < n ; i++)
				for(int j = 0 ; j < n ; j++)
					genC[i][j] = i == j ? 1 : 0;
		}
		
		// Generate new solutions
		ArrayList<SolutionRoot> newRoots = new ArrayList<SolutionRoot>();
		for(int i = 0 ; i < number ; i++)
		{
			int var = 0;
			double[] values = MultiVarNormal.sample(mean, genC);
			ArrayList<Integer> discVals = new ArrayList<Integer>();
			for(int j = 0 ; j < discVars.size() ; j++)
			{
				int value = (int)values[var];
				int min = discVars.get(j).getMin();
				int max = min + discVars.get(j).getCount() - 1;
				value = value < min ? min : value;
				value = value > max ? max : value;
				discVals.add(value);
				var++;
			}
			ArrayList<Double> contVals = new ArrayList<Double>();
			for(int j = 0 ; j < contVars.size() ; j++)
			{
				double value = values[var];
				double min = contVars.get(j).getMin();
				double max = contVars.get(j).getMax();
				value = value < min ? min : value;
				value = value > max ? max : value;
				contVals.add(value);
				var++;
			}
			SolutionRoot root = new SolutionRoot(discVals, contVals);
			newRoots.add(root);
		}
		return newRoots;
	}
	
	/**
	 * Initializes the attributes when a new population size is used
	 * @param n The number of variables
	 * @param popSize The number of solutions in the current population
	 */
	private void initAttributes(int n, int popSize)
	{
		// Initialize attributes
		iteration = 0;
		miu = Math.max(2, (int)((1 - trunc)*popSize));
		computeWeights(miu);
		double sum = 0;
		for(int i = 0 ; i < weights.length ; i++)
			sum += weights[i]*weights[i];
		miuEff = 1/sum;
		cSigma = (miuEff + 2)/(n + miuEff + 3);
		dSigma = 1 + 2*Math.max(0, Math.sqrt((miuEff - 1)/(n + 1)) - 1) + cSigma;
		cC = 4.0/(n + 4.0);
		miuCov = miuEff;
		cCov = 2/(miuCov*(n + 1.4142135623731)*(n + 1.4142135623731)) 
						+ (1 - (1/miuCov))*Math.min(1, (2*miuEff - 1)/((n + 2)*(n + 2) + miuEff));
		
		// Initialize mean and overall standard deviation
		mean = new double[n];
		int var = 0;
		ContSeries stDevs = new ContSeries();
		for(int i = 0 ; i < discVars.size() ; i++)
		{
			double min = discVars.get(i).getMin();
			double max = min + discVars.get(i).getCount() - 1;
			double meanValue = (min + max)/2;
			mean[var] = meanValue;
			stDevs.addValue(Uniform.computeStDev(min, max));
			var++;
		}
		for(int i = 0 ; i < contVars.size() ; i++)
		{
			double min = contVars.get(i).getMin();
			double max = contVars.get(i).getMax();
			double meanValue = (max + min)/2;
			mean[var] = meanValue;
			stDevs.addValue(Uniform.computeStDev(min, max));
			var++;
		}
		sigma = stDevs.getMean();
		
		// Initialize vectors
		pSigma = new double[n];
		pC = new double[n];
		
		// Initialize covariance matrix as the identity matrix
		C = new double[n][n];
		for(int i = 0 ; i < n ; i++)
			for(int j = 0 ; j < n ; j++)
				C[i][j] = i == j ? 1 : 0;
		cSqrtInvLife = Integer.MAX_VALUE;
	}
	
	/**
	 * Computes the weighting factors according to the value of <i>q</i> and the size of the 
	 * population. These weights are used for computing the contribution of each solution in the
	 * construction of the multivariate normal distribution to generate new solutions from.
	 * @param popSize The size of the population
	 */
	private void computeWeights(int popSize)
	{
		weights = null;
		
		// Compute weights
		if(Double.isNaN(q))
		{
			weights = new double[popSize];
			for(int i = 0 ; i < popSize ; i++)
				weights[i] = 1.0/popSize;
		}
		else if(q > 0)
		{
			weights = new double[popSize];
			double sum = 0;
			for(int i = 0 ; i < popSize ; i++)
			{
				int rank = i + 1;
				double weight = 1/(q*popSize*Math.sqrt(2*Math.PI))
								* Math.exp(-(rank - 1)*(rank - 1)/
										(2*q*q*popSize*popSize));
				weights[i] = weight;
				sum += weight;
			}
			
			// Normalize
			for(int i = 0 ; i < popSize ; i++)
				weights[i] = weights[i]/sum;
		}
		
	}
	
}
