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

package maestro.gen.gd;

import java.util.ArrayList;

import maestro.gen.gd.NormSolution;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import probDist.multiVar.MultiVarNormal;
import utilities.stat.ContSeries;

/**
 * Represents a vicinity or neighborhood of solutions for a multi-objective optimization problem. 
 * The vicinity is to be used to find new candidate solutions in the direction of the most 
 * promising measured gradients within a Gradient Descent framework.
 * @author Felipe Hernández
 */
public class Vicinity
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Precision to check for the equality of two double numbers 
	 */
	private static final double DOUBLE_PRECISION = 1E-16;
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The solution on which the vicinity is centered
	 */
	private NormSolution base;
	
	/**
	 * The solutions in the proximity of the {@link #base}. There should be as many neighbors as 
	 * the number of decision variables (<i>n</i>) so that the gradients can be computed. Also, 
	 * they should be able to define, together with the {@link #base}, a basis for the 
	 * <i>n</i>-dimensional space. That is, any three solutions should be linearly independent.
	 */
	private ArrayList<NormSolution> neighbors;
	
	/**
	 * A matrix of size <i>n</i> x <i>n</i> (<i>n</i> is the number of decision variables) with the 
	 * differences between the {@link #neighbors} and the {@link #base}. Rows correspond to 
	 * neighbors and columns correspond to dimensions.
	 */
	private double[][] deltaMat;
	
	/**
	 * A vector with the fitness gradient computed from the solutions in the vicinity
	 */
	private double[] gradient;
	
	// --------------------------------------------------------------------------------------------
	// Constructors
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new vicinity
	 * @param base {@link #base}
	 */
	public Vicinity(NormSolution base)
	{
		this.base		= base;
		int n			= base.getValues().length;
		neighbors		= new ArrayList<>();
		deltaMat		= new double[n][n];
		gradient		= null;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Adds a new solution if it has the adequate number of {@link NormSolution#values} (<i>n</i>) 
	 * and if it contributes to create a basis for the <i>n</i>-dimensional decision variable 
	 * space. That is, any three solutions should be linearly independent. Also, only a number of 
	 * solutions equal to the number of decision variables are accepted.
	 * @param solution The solution to be added to the vicinity
	 * @return True if the solution was added
	 */
	public boolean offerSolution(NormSolution solution)
	{
		int n = base.getValues().length;
		
		// Verify if the neighbors are complete
		if (n < neighbors.size())
			return false;
		
		// Verify number of decision variables
		if (n > solution.getValues().length)
			return false;
		
		// Compute delta vector
		double[] delta0	= new double[n];
		for (int k = 0; k < n; k++)
			delta0[k]	= solution.getValues()[k] - base.getValues()[k];
		
		// Verify the solution is not collinear with any other two solutions
		for (int j = 0; j < neighbors.size(); j++)
		{
			RealVector v1 = new ArrayRealVector(delta0);
			RealVector v2 = new ArrayRealVector(deltaMat[j]);
			if (1.0 - Math.abs(v1.cosine(v2)) < DOUBLE_PRECISION)
				return false;
		}
		for (int i = 0; i < neighbors.size() - 1; i++)
		{
			for (int j = i + 1; j < neighbors.size(); j++)
			{
				NormSolution sol1	= neighbors.get(i);
				NormSolution sol2	= neighbors.get(j);
				double[] delta1		= new double[n];
				double[] delta2		= new double[n];
				for (int k = 0; k < n; k++)
				{
					delta1[k]		= sol2.getValues()[k]		- sol1.getValues()[k];
					delta2[k]		= solution.getValues()[k]	- sol1.getValues()[k];
				}
				RealVector v1		= new ArrayRealVector(delta1);
				RealVector v2		= new ArrayRealVector(delta2);
				if (1.0 - Math.abs(v1.cosine(v2)) < DOUBLE_PRECISION)
					return false;
			}
		}
		
		// Add solution
		neighbors.add(solution);
		deltaMat[neighbors.size() - 1] = delta0;
		return true;
	}
	
	/**
	 * @return The number of solutions in the vicinity including the base
	 */
	public int size()
	{
		return neighbors.size() + 1;
	}
	
	/**
	 * @return {@link #gradMagnitude}
	 */
	public double getGradMagnitude()
	{
		if (gradient == null)
			computeGradient();
		if (gradient == null)
			return Double.NaN;
		RealVector gradVec = new ArrayRealVector(gradient);
		return gradVec.getNorm();
	}
	
	/**
	 * @return {@link #gradient}
	 */
	public double[] getGradient()
	{
		return gradient;
	}
	
	/**
	 * Randomly generates a number of solutions in the direction of the gradient. A Gaussian Kernel
	 * is used to sample the new solutions from. The kernel is centered in the direction of the 
	 * gradient.
	 * @param count 	The number of solutions to generate
	 * @param stepSize 	The factor to multiply the gradient vector to produce the displacement 
	 * vector from the location of the {@link #base} to the the center of the generation kernel 
	 * @param amplitude The factor to multiply the standard deviation (for each dimension) of the 
	 * solutions in the vicinity to set to the diagonal covariance matrix of the Gaussian 
	 * generation kernel
	 * @return The list of randomly generated solutions. <code>null</code> if no solutions can be
	 * generated.
	 */
	public ArrayList<NormSolution> randomGenerate(int count, double stepSize, double amplitude)
	{
		// Get dimensions and verify neighbors
		int n						= base.getValues().length;
		if (neighbors.size() < n)
			return null;
		
		// Compute gradients and generation vectors
		if (gradient == null)
			if (!computeGradient())
				return null;
		double[] diagonal			= new double[n];
		double[] extent				= new double[n];
		for (int i = 0; i < n; i++)
		{
			ContSeries statistics	= new ContSeries(false);
			statistics.addValue(base.getValues()[i]);
			for (int j = 0; j < n; j++)
				statistics.addValue(neighbors.get(j).getValues()[i]);
			double stDev			= statistics.getStDev();
			diagonal[i]				= stDev*stDev*amplitude*amplitude;
			extent[i]				= stDev;
		}
		RealVector origin			= new ArrayRealVector(base.getValues());
		
		// Generate candidate solutions
		RealVector gradVec		= new ArrayRealVector(gradient);		
		RealVector extentVec	= new ArrayRealVector(extent);
		double extentD			= gradVec.dotProduct(extentVec);
		RealVector root			= origin.add(gradVec.mapMultiply(-stepSize*Math.abs(extentD)));
		
		//RealVector root			= origin.add(gradVec.mapMultiply(-stepSize));
		
		int baseIndex			= base.getBaseIndex();
		ArrayList<NormSolution> candidates = new ArrayList<>();
		for (int k = 0; k < count; k++)
		{
			if (amplitude > 0.0)
			{
				MultiVarNormal kernel	= new MultiVarNormal(root.toArray(), diagonal);
				candidates.add(new NormSolution(baseIndex, kernel.sample(), Double.NaN));
			}
			else
				candidates.add(new NormSolution(baseIndex, root.toArray(), Double.NaN));
		}
		return candidates;
	}
	
	/**
	 * Computes {@link #gradient} if the number of solutions in {@link #neighbors} is sufficient
	 * @return True if the gradient was computed. False otherwise (the number of neighbors is
	 * insufficient or the delta matrix cannot be inverted).
	 */
	private boolean computeGradient()
	{
		// Get dimensions and verify neighbors
		int n						= base.getValues().length;
		if (neighbors.size() < n)
			return false;
		
		// Compute gradient
		double[] deltaf			= new double[n];
		for (int j = 0; j < n; j++)
			deltaf[j]			= neighbors.get(j).getFitness() - base.getFitness();
		RealVector deltaFV		= new ArrayRealVector(deltaf);
		try
		{
			RealMatrix deltaRM			= new Array2DRowRealMatrix(deltaMat).transpose();
			DecompositionSolver solver	= new LUDecomposition(deltaRM).getSolver();
			RealVector gradVec		= solver.solve(deltaFV);
			gradient				= gradVec.toArray();
		} catch (Exception e)
		{
			gradient				= null;
			return false;
		}
		return true;
	}
	
}
