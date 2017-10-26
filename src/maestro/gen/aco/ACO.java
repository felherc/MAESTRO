package maestro.gen.aco;

import java.util.ArrayList;

import maestro.ContVar;
import maestro.DiscVar;
import maestro.gen.Generator;
import maestro.gen.aco.ContinuousVariables;
import maestro.gen.aco.DiscreteVariables;
import maestro.solution.Solution;
import maestro.solution.SolutionRoot;

/**
 * This generator allows to create new solutions to a problem with discrete and continuous 
 * variables using a sorted base population of solutions and an Ant Colony Optimization algorithm
 * @author Felipe Hernández
 */
public class ACO implements Generator 
{
	
	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The identifier of the generator method
	 */
	public final static String ID = "Ant Colony Optimization";
	
	/**
	 * The short version of the identifier of the generator method
	 */
	public final static String SHORT_ID = "ACO";
	
	/**
	 * The <code>minPheromone</code> attribute printing name
	 */
	public final static String PARAM_MIN_PHEROMONE = "Minimum pheromone = ";
	
	/**
	 * The <code>elitism</code> attribute printing name
	 */
	public final static String PARAM_ELITISM = "Elitism = ";
	
	/**
	 * The <code>randomP</code> attribute printing name
	 */
	public final static String PARAM_RANDOM_P = "Random probability = ";
	
	/**
	 * The <code>q</code> attribute printing name
	 */
	public final static String PARAM_Q = "q = ";
	
	/**
	 * The <code>xi</code> attribute printing name
	 */
	public final static String PARAM_XI = "xi = ";
	
	/**
	 * Default value for the minimum amount of pheromone a value can have expressed as a percentage 
	 * of the average pheromone value.
	 */
	public final static double D_MIN_PHEROMONE = 0.3; 
	// Does not make much difference
	
	/**
	 * Default value for the rate at which the best solution found so far deposits pheromone. 
	 * Preferably between 0 and 1.
	 */
	public final static double D_ELITISM = 0.0; 
	// Makes small difference
	
	/**
	 * Default value for the probability of generating uniformly-distributed random values for 
	 * continuous variables in new solutions
	 */
	public final static double D_RANDOM_P = 0.1; 
	// Very important - sweet spot around 0.05-0.15; smaller is better for higher dimensionality
	
	/**
	 * Default value for the weight assignment variable q
	 */
	public final static double D_Q = 0.3; 
	// Relatively indifferent for range 0.2-0.8 for few dimensions; make it smaller for high dimensionality
	
	/**
	 * Default value for the standard deviation computing parameter xi
	 */
	public final static double D_XI = 0.5; 
	// Very important - sweet-spot around 0.4-0.5; higher for higher dimensionality

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Object that manages the discrete variables of the problem
	 */
	private DiscreteVariables discreteVariables;
	
	/**
	 * Object that manages the continuous variables of the problem
	 */
	private ContinuousVariables continuousVariables;
	
	/**
	 * The minimum amount of pheromone a value can have expressed as a percentage of the average 
	 * pheromone value 
	 */
	private double minPheromone;
	
	/**
	 * The rate at which the best solution found so far deposits pheromone. The complement is 
	 * deposited by the best ants in the iteration. Preferably between 0 and 1.
	 */
	private double elitism;
	
	/**
	 * The probability of generating uniformly-distributed random values for 
	 * continuous variables in new solutions
	 */
	private double randomP;
	
	/**
	 * This variable defines how the solutions in the buffer are assigned weights. The 
	 * weights are used to compute the pheromone update. Q has to be positive. A value of 0 means 
	 * only the most fitted solution may lay pheromone. A large value of q means that every 
	 * solution contributes to the pheromone with little regard to its ranking.
	 */
	private double q;
	
	/**
	 * This parameter is used to compute the standard deviation of normal probability 
	 * functions in order to generate new values for the variables. When generating a new value for 
	 * a given variable, a Gaussian kernel is created from the values of the solutions in the 
	 * archive. A normal distribution is created for each solution with its value being the mean 
	 * and an standard deviation. This deviation is computed by multiplying the xi parameter with 
	 * the average distance between the mean and the values of the other solutions (the other 
	 * means). xi must be positive. A small value tends to speed convergence, while a large value 
	 * tends to perform a more complete search of the solution space.
	 */
	private double xi;
	
	// --------------------------------------------------------------------------------------------
	// Constructors
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new ACO generator instance with the default parameter values
	 */
	public ACO()
	{
		discreteVariables = null;
		continuousVariables = null;
		minPheromone = D_MIN_PHEROMONE;
		elitism = D_ELITISM;
		randomP = D_RANDOM_P;
		q = D_Q;
		xi = D_XI;
	}
	
	/**
	 * Creates a new ACO generator instance
	 * @param minPheromone The minimum amount of pheromone a value can have expressed as a 
	 * percentage of the average pheromone value
	 * @param elitism The rate at which the best solution found so far deposits pheromone. The 
	 * complement is deposited by the best ants in the iteration. Preferably between 0 and 1.
	 * @param randomP The probability of generating uniformly-distributed random values for 
	 * continuous variables in new solutions
	 * @param q This variable defines how the solutions in the buffer are assigned weights. The 
	 * weights are used to compute the pheromone update. Q has to be positive. A value of 0 means 
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
	public ACO(double minPheromone, double elitism, double randomP, double q, double xi)
	{
		discreteVariables = null;
		continuousVariables = null;
		this.minPheromone = minPheromone;
		this.elitism = elitism;
		this.randomP = randomP;
		this.q = q;
		this.xi = xi;
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
	 * @return The minimum amount of pheromone a value can have expressed as a percentage of the average 
	 * pheromone value
	 */
	public double getMinPheromone() 
	{
		return minPheromone;
	}

	/**
	 * @param minPheromone the minPheromone to set
	 */
	public void setMinPheromone(double minPheromone) 
	{
		this.minPheromone = minPheromone;
	}

	/**
	 * @return The rate at which the best solution found so far deposits pheromone. The complement 
	 * is deposited by the best ants in the iteration. Preferably between 0 and 1.
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
		this.elitism = elitism > 1 ? 1 : elitism;
	}

	/**
	 * @return The probability of generating uniformly-distributed random values for continuous 
	 * variables in new solutions
	 */
	public double getRandomP() 
	{
		return randomP;
	}

	/**
	 * @param randomP The probability of generating uniformly-distributed random values for 
	 * continuous variables in new solutions
	 */
	public void setRandomP(double randomP) 
	{
		this.randomP = randomP < 0 ? 0 : randomP;
		this.randomP = randomP > 1 ? 1 : randomP;
	}

	/**
	 * @return This variable defines how the solutions in the buffer are assigned weights. The 
	 * weights are used to compute the pheromone update. Q has to be positive. A value of 0 means 
	 * only the most fitted solution may lay pheromone. A large value of q means that every 
	 * solution contributes to the pheromone with little regard to its ranking.
	 */
	public double getQ() 
	{
		return q;
	}

	/**
	 * @param q This variable defines how the solutions in the buffer are assigned weights. The 
	 * weights are used to compute the pheromone update. Q has to be positive. A value of 0 means 
	 * only the most fitted solution may lay pheromone. A large value of q means that every 
	 * solution contributes to the pheromone with little regard to its ranking.
	 */
	public void setQ(double q) 
	{
		this.q = q > 0 ? 0 : q;
	}

	/**
	 * @return This parameter is used to compute the standard deviation of normal probability 
	 * functions in order to generate new values for the variables. When generating a new value for 
	 * a given variable, a Gaussian kernel is created from the values of the solutions in the 
	 * archive. A normal distribution is created for each solution with its value being the mean 
	 * and an standard deviation. This deviation is computed by multiplying the xi parameter with 
	 * the average distance between the mean and the values of the other solutions (the other 
	 * means). xi must be positive. A small value tends to speed convergence, while a large value 
	 * tends to perform a more complete search of the solution space.
	 */
	public double getXi() 
	{
		return xi;
	}

	/**
	 * @param xi This parameter is used to compute the standard deviation of normal probability 
	 * functions in order to generate new values for the variables. When generating a new value for 
	 * a given variable, a Gaussian kernel is created from the values of the solutions in the 
	 * archive. A normal distribution is created for each solution with its value being the mean 
	 * and an standard deviation. This deviation is computed by multiplying the xi parameter with 
	 * the average distance between the mean and the values of the other solutions (the other 
	 * means). xi must be positive. A small value tends to speed convergence, while a large value 
	 * tends to perform a more complete search of the solution space.
	 */
	public void setXi(double xi) 
	{
		this.xi = xi < 0 ? 0 : xi;
	}
	
	@Override
	public String getParamSummary() 
	{
		String line = "";
		line += PARAM_MIN_PHEROMONE + minPheromone;
		line += "; " + PARAM_ELITISM + elitism;
		line += "; " + PARAM_RANDOM_P + randomP;
		line += "; " + PARAM_Q + q;
		line += "; " + PARAM_XI + xi;
		return line;
	}

	@Override
	public void addDiscVariable(DiscVar variable) 
	{
		if(discreteVariables == null)
			discreteVariables = new DiscreteVariables(minPheromone, elitism, q);
		String name = variable.getName();
		ArrayList<String> values = variable.getValues();
		if(values == null)
		{
			int min = variable.getMin();
			int valueCount = variable.getCount();
			discreteVariables.addVariable(name, min, valueCount);
		}
		else
			discreteVariables.addVariable(name, values);
	}

	@Override
	public void addContVariable(ContVar variable) 
	{
		String name = variable.getName();
		double min = variable.getMin();
		double max = variable.getMax();
		if(continuousVariables == null)
			continuousVariables = new ContinuousVariables(1 - randomP, q, xi);
		continuousVariables.addVariable(name, min, max);
	}

	@Override
	public void clearVariables() 
	{
		discreteVariables = null;
		continuousVariables = null;
	}

	@Override
	public synchronized ArrayList<SolutionRoot> generateSolutions(
													ArrayList<Solution> population, int number)
	{	
		if(population.size() > 0 || discreteVariables != null || continuousVariables != null)
		{
			// Reset values and adjust parameters
			if(discreteVariables != null)
			{
				discreteVariables.reset();
				discreteVariables.setMinPheromone(minPheromone);
				discreteVariables.setElitism(elitism);
				discreteVariables.setQ(q);
			}
			if(continuousVariables != null)
			{
				continuousVariables.reset();
				continuousVariables.setElitism(1 - randomP);
				continuousVariables.setQ(q);
				continuousVariables.setXi(xi);
			}
			
			// Set pheromone levels
			for(Solution solution : population)
			{
				if(discreteVariables != null)
					discreteVariables.addSolution(solution);
				if(continuousVariables != null)
					continuousVariables.addSolution(solution);
			}
			
			// Create solutions
			ArrayList<SolutionRoot> newSolutions = new ArrayList<SolutionRoot>();
			ArrayList<Integer> discValues = null;
			ArrayList<Double> contValues = null;
			for(int i = 0 ; i < number ; i++)
			{
				if(discreteVariables != null)
					discValues = discreteVariables.generateSolution();
				if(continuousVariables != null)
					contValues = continuousVariables.generateSolution();
				SolutionRoot newSolution = new SolutionRoot(discValues, contValues);
				newSolutions.add(newSolution);
			}
			return newSolutions;
		}
		return null;
	}

}
