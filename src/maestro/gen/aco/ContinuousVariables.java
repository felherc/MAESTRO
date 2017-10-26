package maestro.gen.aco;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;

import maestro.gen.aco.ContinuousVariable;
import maestro.solution.Solution;
import probDist.Normal;

/**
 * This class stores the information of the continuous variables of an Ant Colony Optimization 
 * problem. This information includes the name of the variables and the pheromone level for each. 
 * New continuous values for a generated solution can be generated according to pheromone values, 
 * and these can be updated by introducing new candidate solutions.
 * @author Felipe Hernández
 */
public class ContinuousVariables
{
	
	// -------------------------------------------------------------------------------------------- 
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * This tree keeps the best solutions found so far in ascending fitness order
	 */
	private TreeSet<Solution> archive;
	
	/**
	 * The list with the names and limit values for the continuous variables
	 */
	private ArrayList<ContinuousVariable> variables;
	
	/**
	 * Stores the weighting factors for computing pheromone contribution of the solutions in the 
	 * archive according to their rank. For maximization problems, weights are assigned in 
	 * ascending order. Conversely, in minimization problems, weights are assigned in descending 
	 * order. Weights follow a Gaussian function using the rank as the parameter, a mean of 1.0 and 
	 * a deviation of k*q, where k is the number of solutions and q is the weighting attribute.
	 */
	private ArrayList<Double> weights;
	
	/**
	 * The probability of selecting a solution in the archive to generate values for new solutions. 
	 * The complementary probability corresponds to generating uniformly distributed values. The 
	 * selection of the generation method is done once per variable per solution generated.
	 */
	private double elitism;
	
	/**
	 * This parameter defines how the solutions in the archive are assigned weights. The weights 
	 * are used to compute the pheromone update. Q has to be positive. A value of 0 means only the 
	 * most fitted solution may lay pheromone. A large value of q means that every solution 
	 * contributes to the pheromone with little regard for their fitness value.
	 */
	private double q;
	
	/**
	 * This parameter is used to compute the standard deviation of normal probability functions in 
	 * order to generate new values for the variables. When generating a new value for a given 
	 * variable, a Gaussian kernel is created from the values of the solutions in the archive. A 
	 * normal distribution is created for each solution with its value being the mean and an 
	 * standard deviation. This deviation is computed by multiplying the xi parameter with the 
	 * average distance between the mean and the values of the other solutions (the other means). 
	 * xi must be positive. A small value tends to speed convergence, while a large value tends to 
	 * perform a more complete search of the solution space.
	 */
	private double xi;
	
	// -------------------------------------------------------------------------------------------- 
	// Constructors
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @param elitism The probability of selecting a solution in the archive to generate values for 
	 * new solutions. The complementary probability corresponds to generating uniformly distributed 
	 * values. The selection of the generation method is done once per variable per solution 
	 * generated. Preferably between 0 and 1.
	 * @param q This parameter defines how the solutions in the archive are assigned weights. The 
	 * weights are used to compute the pheromone update. q has to be positive. A value of 0 means 
	 * only the most fitted solution may lay pheromone. A large value of q means that every 
	 * solution contributes to the pheromone with little regard for their fitness value.
	 * @param xi This parameter is used to compute the standard deviation of normal probability 
	 * functions in order to generate new values for the variables. When generating a new value for 
	 * a given variable, a Gaussian kernel is created from the values of the solutions in the 
	 * archive. A normal distribution is created for each solution with its value being the mean 
	 * and an standard deviation. This deviation is computed by multiplying the xi parameter with 
	 * the average distance between the mean and the values of the other solutions (the other 
	 * means). xi must be positive. A small value tends to speed convergence, while a large value 
	 * tends to perform a more complete search of the solution space.
	 */
	public ContinuousVariables(double elitism, double q, double xi)
	{
		archive = new TreeSet<Solution>();
		variables = new ArrayList<ContinuousVariable>();
		this.elitism = elitism < 0 ? 0 : elitism;
		this.elitism = elitism > 1 ? 1 : this.elitism;
		this.q = q < 0 ? 0 : q;
		this.xi = xi < 0 ? 0 : xi;
		computeWeights();
	}
	
	// -------------------------------------------------------------------------------------------- 
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return The probability of selecting a solution in the archive to generate values for new 
	 * solutions. The complementary probability corresponds to generating uniformly distributed 
	 * values. The selection of the generation method is done once per variable per solution 
	 * generated.
	 */
	public double getElitism() 
	{
		return elitism;
	}

	/**
	 * @param elitism The probability of selecting a solution in the archive to generate values 
	 * for new solutions. The complementary probability corresponds to generating uniformly 
	 * distributed values. The selection of the generation method is done once per variable per 
	 * solution generated. Preferably between 0 and 1.
	 */
	public void setElitism(double elitism) 
	{
		this.elitism = elitism < 0 ? 0 : elitism;
		this.elitism = elitism > 1 ? 1 : this.elitism;
	}

	/**
	 * @return How the solutions in the archive are assigned weights. The weights are used to 
	 * compute the pheromone update. A value of 0 means only the most fitted solution may lay 
	 * pheromone. A large value of q means that every solution contributes to the pheromone with 
	 * little regard for their fitness value.
	 */
	public double getQ()
	{
		return q;
	}

	/**
	 * @param q This parameter defines how the solutions in the archive are assigned weights. The 
	 * weights are used to compute the pheromone update. Q has to be positive. A value of 0 means 
	 * only the most fitted solution may lay pheromone. A large value of q means that every 
	 * solution contributes to the pheromone with little regard for their fitness value.
	 */
	public void setQ(double q) 
	{
		this.q = q < 0 ? 0 : q;
	}

	/**
	 * @return The parameter used to compute the standard deviation of normal probability functions 
	 * in order to generate new values for the variables. When generating a new value for a given 
	 * variable, a Gauss Kernel is created from the values of the solutions in the archive. A 
	 * normal distribution is created for each solution with its value being the mean and an 
	 * standard deviation. This deviation is computed by multiplying the xi parameter with the 
	 * average distance between the mean and the values of the other solutions (the other means). A 
	 * small value tends to speed convergence, while a large value tends to perform a more complete 
	 * search of the solution space.
	 */
	public double getXi() 
	{
		return xi;
	}

	/**
	 * @param xi This parameter is used to compute the standard deviation of normal probability 
	 * functions in order to generate new values for the variables. When generating a new value for 
	 * a given variable, a Gauss Kernel is created from the values of the solutions in the archive. 
	 * A normal distribution is created for each solution with its value being the mean and an 
	 * standard deviation. This deviation is computed by multiplying the xi parameter with the 
	 * average distance between the mean and the values of the other solutions (the other means). 
	 * xi must be positive. A small value tends to speed convergence, while a large value tends to 
	 * perform a more complete search of the solution space.
	 */
	public void setXi(double xi) 
	{
		this.xi = xi < 0 ? 0 : xi;
	}
	
	/**
	 * Creates a new continuous variable
	 * @param name The identifier of the variable
	 * @param minValue The minimum value the variable can take
	 * @param maxValue The maximum value the variable can take
	 * @return The index of the created variable
	 * @throws Exception If the maximum value is not greater than the minimum value
	 */
	public int addVariable(String name, double minValue, double maxValue)
	{
		ContinuousVariable variable = new ContinuousVariable(name, minValue, maxValue);
		variables.add(variable);
		return variables.size() - 1;
	}
	
	/**
	 * Adds a new generated solution to the archive. If the archive reaches the maximum, the best 
	 * fitted solution is removed. This way, the archive stores the best solutions that are then 
	 * used to generate new values.
	 * @param solution The solution to add to the archive
	 */
	public void addSolution(Solution solution)
	{
		if(solution.getContValues().size() == variables.size())
			archive.add(solution);
	}
	
	/**
	 * @return A list with the double values of a generated solution. The values are generated by 
	 * sampling a Gaussian kernel constructed with the values of the solutions in the archive, the 
	 * weighting factors and the xi parameter. If the archive has not reached its maximum size, 
	 * uniformly distributed values within each variable's range are returned.
	 */
	public synchronized ArrayList<Double> generateSolution()
	{
		computeWeights();
		ArrayList<Double> generatedValues = new ArrayList<Double>();
		if(archive.size() == 0)
			return generatedValues;
		
		for(int i = 0 ; i < variables.size() ; i++)
		{
			double value = 0;
			double elitist = Math.random();
			if(elitist >= elitism)
			{
				ContinuousVariable variable = variables.get(i);
				double min = variable.getMinValue();
				double max = variable.getMaxValue();
				value = min + Math.random()*(max - min);
			}
			else
			{
				// Obtain solution index
				int index = -1;
				double sum = 0;
				if(weights != null)
				{
					double random = Math.random();
					for(int j = 0 ; j < weights.size() ; j++)
					{
						sum += weights.get(j);
						if(sum >= random)
						{
							index = j;
							j = weights.size(); // Breaks the cycle
						}
					}
				}
				else
					index = archive.size() - 1;
				
				// Get mean value
				Iterator<Solution> iterator = archive.iterator();
				int j = 0;
				Solution solution = iterator.next();
				while(j < index)
				{
					solution = iterator.next();
					j++;
				}
				double mean = solution.getContValues().get(i);
				
				// Compute standard deviation
				iterator = archive.iterator();
				sum = 0;
				for(j = 0 ; j < archive.size() ; j++)
				{
					solution = iterator.next();
					value = solution.getContValues().get(i);
					sum += Math.abs(mean - value);
				}
				double deviation = xi*sum/(archive.size() - 1);
				
				// Generate random value and add to solution list
				boolean ok = false;
				while(!ok)
				{
					value = Normal.sample(mean, deviation);
					ContinuousVariable variable = variables.get(i);
					if(value >= variable.getMinValue() && 
										value <= variable.getMaxValue())
						ok = true;
				}
			}
			generatedValues.add(value);
		}
		
		return generatedValues;
	}
	
	/**
	 * Initializes the weighting factors according to the value of q and the maximum size of the 
	 * archive. These weights are used for computing the probability of selecting a given solution 
	 * in the archive to generate new values for the variables.
	 */
	private void computeWeights()
	{
		if(q > 0)
		{
			// Compute weights
			weights = new ArrayList<Double>();
			double sum = 0;
			int archiveMaxSize = archive.size();
			for(int i = 0 ; i < archiveMaxSize ; i++)
			{
				int rank = archiveMaxSize - i;
				double weight = 1/(q*archiveMaxSize*Math.sqrt(2*Math.PI))
								* Math.exp(-(rank - 1)*(rank - 1)/
										(2*q*q*archiveMaxSize*archiveMaxSize));
				weights.add(weight);
				sum += weight;
			}
			
			// Normalize
			for(int i = 0 ; i < archiveMaxSize ; i++)
				weights.set(i, weights.get(i)/sum);
		}
		else
			weights = null;
	}
	
	/**
	 * Returns the name of the variable
	 * @param index The index of the variable
	 * @return The name of the variable
	 */
	public String getName(int index)
	{
		return variables.get(index).getName();
	}
	
	/**
	 * Clears the archive with the best solutions found so far
	 */
	public void reset()
	{
		archive.clear();
	}
	
}
