package maestro.gen.aco;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;

import maestro.gen.aco.DiscreteVariable;
import maestro.solution.Solution;

/**
 * This class stores the information of the discrete variables of an Ant Colony Optimization 
 * problem. This information includes the name of the variables and the pheromone level for each. 
 * New discrete values for a generated solution can be generated according to pheromone values, and 
 * these can be updated by introducing new candidate solutions.
 * @author Felipe Hernández
 */
public class DiscreteVariables
{
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The list with the discrete variables objects. Each variable includes the pheromone levels to 
	 * generate new solutions.
	 */
	private ArrayList<DiscreteVariable> variables;
	
	/**
	 * The best solution found so far
	 */
	private Solution best;
	
	/**
	 * The rate at which the best solution found so far deposits pheromone. The complement is 
	 * deposited by the best ants in the iteration. Preferably between 0 and 1.
	 */
	private double elitism;
	
	/**
	 * Ordered structure with the most recently created ant solutions. The solutions are ordered 
	 * from worst to best according to their fitness value. New solutions are added to the buffer 
	 * until it reaches the capacity specified by the bufferSize attribute. Then, the pheromone 
	 * levels are updated and the buffer is emptied.
	 */
	private TreeSet<Solution> buffer;
	
	/**
	 * Stores the weighting factors for computing pheromone contribution of the solutions in the 
	 * buffer according to their rank. For maximization problems, weights are assigned in ascending 
	 * order. Conversely, in minimization problems, weights are assigned in descending order. 
	 * Weights follow a Gaussian function using the rank as the parameter, a mean of 1.0 and a 
	 * deviation of k*q, where k is the number of solutions and q is the weighting attribute. If 
	 * the instance is equal to null, only the best ranked solution deposits pheromone.
	 */
	private ArrayList<Double> weights;
	
	/**
	 * This parameter defines how the solutions in the buffer are assigned weights. The weights are 
	 * used to compute the pheromone update. Q has to be positive. A value of 0 means only the most 
	 * fitted solution may lay pheromone. A large value of q means that every solution contributes 
	 * to the pheromone with little regard for their fitness value.
	 */
	private double q;
	
	/**
	 * The minimum amount of pheromone a value can have expressed as a percentage of the average 
	 * pheromone value. Used to initialize new variables. 
	 */
	private double minPheromone;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @param minPheromone The minimum amount of pheromone a value can have expressed as a 
	 * percentage of the average pheromone value. Used to initialize new variables.
	 * @param elitism The rate at which the best solution found so far deposits pheromone. The 
	 * complement is deposited by the best ants in the iteration. Preferably between 0 and 1.
	 * @param q This variable defines how the solutions in the buffer are assigned weights. The 
	 * weights are used to compute the pheromone update. Q has to be positive. A value of 0 means 
	 * only the most fitted solution may lay pheromone. A large value of q means that every 
	 * solution contributes to the pheromone with little regard for their fitness value.
	 */
	public DiscreteVariables(double minPheromone, double elitism, double q)
	{
		variables = new ArrayList<DiscreteVariable>();
		this.minPheromone = minPheromone;
		best = null;
		this.elitism = elitism < 0 ? 0 : elitism;
		this.elitism = elitism > 1 ? 1 : this.elitism;
		buffer = new TreeSet<Solution>();
		this.q = q < 0 ? 0 : q;
		computeWeights();
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	/**
	 * @return The rate at which the best solution found so far deposits pheromone. The complement 
	 * is deposited by the best ants in the iteration.
	 */
	public double getElitism() 
	{
		return elitism;
	}

	/**
	 * @param elitism The rate at which the best solution found so far deposits pheromone. The 
	 * complement is deposited by the best ants in the iteration. Preferably between 0 and 1.
	 */
	public void setElitism(double elitism) 
	{
		this.elitism = elitism < 0 ? 0 : elitism;
		this.elitism = elitism > 1 ? 1 : this.elitism;
	}

	/**
	 * @return How the solutions in the buffer are assigned weights. Q has to be positive. A value 
	 * of 0 means only the most fitted solution may lay pheromone. A large value of q means that 
	 * every solution contributes to the pheromone with little regard for their fitness value.
	 */
	public double getQ()
	{
		return q;
	}

	/**
	 * @param q This variable defines how the solutions in the buffer are assigned weights. The 
	 * weights are used to compute the pheromone update. Q has to be positive. A value of 0 means 
	 * only the most fitted solution may lay pheromone. A large value of q means that every 
	 * solution contributes to the pheromone with little regard for their fitness value.
	 */
	public void setQ(double q)
	{
		this.q = q < 0 ? 0 : q;
	}
	
	/**
	 * @param percentage The percentage of the average pheromone per value to be set as the minimum 
	 * pheromone value. Preferably between 0 and 1.
	 */
	public void setMinPheromone(double percentage)
	{
		for(DiscreteVariable variable : variables)
			variable.setMinPheromone(percentage);
	}
	
	/**
	 * Creates a new discrete variable
	 * @param name The identifier of the variable
	 * @param min The minimum integer value the variable can take
	 * @param valueCount The number of values the variable can take
	 * @return The index of the created variable
	 */
	public int addVariable(String name, int min, int valueCount)
	{
		DiscreteVariable variable = new DiscreteVariable(name, min, valueCount);
		variable.setMinPheromone(minPheromone);
		variables.add(variable);
		return variables.size() - 1;
	}
	
	/**
	 * Creates a new discrete variable
	 * @param name The identifier of the variable
	 * @param values The set of identifiers of the variable's possible values
	 * @return The index of the created variable
	 */
	public int addVariable(String name, ArrayList<String> values)
	{
		DiscreteVariable variable = new DiscreteVariable(name, values);
		variable.setMinPheromone(minPheromone);
		variables.add(variable);
		return variables.size() - 1;
	}
	
	/**
	 * Adds a new generated solution to the buffer. If the buffer reaches the maximum size, the 
	 * pheromone levels are updated and the buffer is cleared. If the solution has a better fitness 
	 * value than the best so far, the solution becomes the new best.
	 * @param solution The solution to add to the buffer
	 */
	public void addSolution(Solution solution)
	{
		// Update best and add to buffer
		if(solution.getDiscValues().size() == variables.size())
		{
			if(best == null)
				best = solution;
			else
			{
				if(solution.compareTo(best) > 0)
					best = solution;
			}
			buffer.add(solution);
		}
	}
	
	/**
	 * @return A list with the integer values of a generated solution. The generation is done 
	 * according to the pheromone levels of each variable.
	 */
	public synchronized ArrayList<Integer> generateSolution()
	{
		computeWeights();
		updatePheromone();
		ArrayList<Integer> solution = new ArrayList<Integer>();
		for(DiscreteVariable variable : variables)
			solution.add(variable.generateValue());
		return solution;
	}
	
	/**
	 * Initializes the weighting factors according to the value of q and the size of the buffer. 
	 * These weights are used for computing pheromone contribution of the solutions in the buffer 
	 * according to their rank.
	 */
	private void computeWeights()
	{
		if(q > 0)
		{
			// Compute weights
			weights = new ArrayList<Double>();
			double sum = 0;
			int bufferSize = buffer.size();
			for(int i = 0 ; i < bufferSize ; i++)
			{
				int rank = bufferSize - i;
				double weight = 1/(q*bufferSize*Math.sqrt(2*Math.PI))
								* Math.exp(-(rank - 1)*(rank - 1)/
										(2*q*q*bufferSize*bufferSize));
				weights.add(weight);
				sum += weight;
			}
			
			// Normalize
			for(int i = 0 ; i < bufferSize ; i++)
				weights.set(i, weights.get(i)/sum);
		}
		else
			weights = null;
	}
	
	/**
	 * Updates the pheromone levels for all the variables using the solutions in the buffer and the 
	 * best solution. The pheromone is updated using the weights established, the pheromone rate, 
	 * the pheromone conservation and the elitism factor. After updating the pheromone values, the 
	 * buffer is cleared.
	 */
	private void updatePheromone()
	{
		for(int i = 0 ; i < variables.size() ; i++)
		{
			// Evaporate pheromone
			DiscreteVariable variable = variables.get(i);
			double pheromone = variable.evaporatePheromone(1.0);
			
			// Compute buffered solutions' pheromone update
			Iterator<Solution> iterator = null;
			if(weights == null)
			{
				iterator = buffer.descendingIterator();
				Solution solution = iterator.next();
				int value = solution.getDiscValues().get(i);
				variable.depositPheromone(value, pheromone*(1 - elitism));
			}
			else
			{
				iterator = buffer.iterator();
				for(int j = 0 ; j < buffer.size() ; j++)
				{
					Solution solution = iterator.next();
					int value = solution.getDiscValues().get(i);
					double deposit = weights.get(j)*pheromone*(1 - elitism);
					variable.depositPheromone(value, deposit);
				}
			}
			
			// Compute best solution's contribution
			if(elitism > 0)
			{
				int value = best.getDiscValues().get(i);
				variable.depositPheromone(value, pheromone*(elitism));
			}
		}
		buffer.clear();
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
	 * Returns the identifier of the specified value of the variable or the value itself if actual 
	 * numbers are being used
	 * @param variableIndex The index of the variable
	 * @param value The index of the value or the integer value of the variable when actual numbers 
	 * are being used
	 * @return The identifier of the specified value of the variable or the value itself if actual 
	 * numbers are being used
	 */
	public String getValueIdentifier(int variableIndex, int value)
	{
		return variables.get(variableIndex).getIdentifier(value);
	}
	
	/**
	 * Resets the pheromone levels to the default values. Clears the best solution found so far and 
	 * the buffer.
	 */
	public void reset()
	{
		for(DiscreteVariable variable : variables)
			variable.reset();
		best = null;
		buffer.clear();
	}
	
}