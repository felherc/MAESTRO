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

import maestro.gen.gd.NormSolution;

/**
 * Represents a normalized solution of a multi-objective optimization problem. Both the values of 
 * the decision variables and the values of the objective functions are normalized to be used in
 * a Gradient Descent framework to find promising exploration gradients to generate new candidate
 * solutions.
 * @author Felipe Hernández
 */
public class NormSolution 
{
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The index of the base of the vicinity in the original solution list. This solution could
	 * either be the base of the vicinity or had been created by it.
	 */
	private int baseIndex;
	
	/**
	 * The normalized values for the decision variables
	 */
	private double[] values;
	
	/**
	 * The normalized value for the objective function. It is assumed that the objective is to be
	 * minimized.
	 */
	private double fitness;
	
	// --------------------------------------------------------------------------------------------
	// Constructors
	// --------------------------------------------------------------------------------------------

	/**
	 * @param values	{@link #values}
	 * @param fitness	{@link #fitness}
	 */
	public NormSolution(double[] values, double fitness) 
	{
		baseIndex		= -1;
		this.values		= values;
		this.fitness	= fitness;
	}
	
	/**
	 * @param vicinityIndex	{@link #baseIndex}
	 * @param values		{@link #values}
	 * @param fitness		{@link #fitness}
	 */
	public NormSolution(int vicinityIndex, double[] values, double fitness) 
	{
		this.baseIndex = vicinityIndex;
		this.values = values;
		this.fitness = fitness;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------	

	/**
	 * @return {@link #baseIndex}
	 */
	public int getBaseIndex() 
	{
		return baseIndex;
	}

	/**
	 * @param vicinityIndex {@link #baseIndex}
	 */
	public void setBaseIndex(int baseIndex) 
	{
		this.baseIndex = baseIndex;
	}

	/**
	 * @return {@link #values}
	 */
	public double[] getValues()
	{
		return values;
	}

	/**
	 * @param values {@link #values}
	 */
	public void setValues(double[] values)
	{
		this.values = values;
	}

	/**
	 * @return {@link #fitness}
	 */
	public double getFitness()
	{
		return fitness;
	}

	/**
	 * @param fitness {@link #fitness}
	 */
	public void setFitness(double fitness)
	{
		this.fitness = fitness;
	}
	
	/**
	 * Computes the Euclidean distance between the two solutions in terms of their {@link #values}
	 * @param other The solution to compute the distance from
	 * @return The Euclidean distance between the two solutions in terms of their {@link #values}
	 */
	public double computeDistance(NormSolution other)
	{
		double sum = 0;
		for (int i = 0; i < values.length; i++)
		{
			double diff		= values[i] - other.getValues()[i];
			sum				+= diff*diff;
		}
		return Math.sqrt(sum);
	}
	
}
