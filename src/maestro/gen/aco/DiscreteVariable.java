package maestro.gen.aco;

import java.util.ArrayList;

/**
 * This class stores the pheromone and other information of a discrete variable to be optimized in 
 * an Ant Colony Optimization problem. Pheromone levels are assigned to each possible value of the 
 * variable and are set to be proportional to the adequacy of that value from the solutions tried 
 * so far. The pheromone levels are used to stochastically generate values for new ants.
 * @author Felipe Hernández
 */
public class DiscreteVariable
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Default value for the pheromone level of variables
	 */
	private final static double DEFAULT_PHEROMONE = 0.5;
	
	/**
	 * Default value for the minimum amount of pheromone a value can have
	 */
	private final static double DEFAULT_MIN_PHEROMONE = 0.1;
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Identifier of the variable
	 */
	private String name;
	
	/**
	 * A list with the pheromone levels for each one of the possible values the variable can take. 
	 * Each level represents a relative measure of the probability of that value being chosen when 
	 * generating new ants. The values must be greater than zero.
	 */
	private ArrayList<Double> pheromone;
	
	/**
	 * The sum of the pheromone levels
	 */
	private double pheromoneSum;
	
	/**
	 * The minimum amount of pheromone a value can have. Must be greater than to zero
	 */
	private double minPheromone;
	
	/**
	 * The minimum value the variable can take. This attribute is used when the variable represents 
	 * an actual number and not just an identifier. The index of the pheromone attribute is added 
	 * to this value to compute the actual values of the variable. Conversely, the attribute values
	 * are used when the variable represents identifiers. If the attribute values is not null, the 
	 * identifier representation is used.
	 */
	private int min;
	
	/**
	 * The identifiers of the values the variable can take. This attribute is used when the values 
	 * do not represent actual numbers but simply identifiers. Use the attribute min to represent 
	 * numbers instead. The instance is null when the latter representation is used.
	 */
	private ArrayList<String> values;
	
	// --------------------------------------------------------------------------------------------
	// Constructors
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Initializes a new discrete variable with its possible values. The values are the specified 
	 * number of integers starting from the min value and sequentially on. Default values for the 
	 * pheromone levels and the minimum pheromone are used.
	 * @param name The identifier of the variable
	 * @param min The minimum integer value the variable can take
	 * @param valueCount The number of values the variable can take
	 */
	public DiscreteVariable(String name, int min, int valueCount)
	{
		this.name = name;
		pheromone = new ArrayList<Double>();
		for(int i = 0 ; i < valueCount ; i++)
			pheromone.add(DEFAULT_PHEROMONE);
		pheromoneSum = DEFAULT_PHEROMONE*pheromone.size();
		minPheromone = DEFAULT_MIN_PHEROMONE;
		this.min = min;
		values = null;
	}
	
	/**
	 * Initializes a new discrete variable with its possible values. The values are integers from 0 
	 * and on, that map to a set of identifiers. Default values for the pheromone levels and the 
	 * minimum pheromone are used. 
	 * @param name The identifier of the variable
	 * @param values The set of identifiers of the variable's possible values
	 */
	public DiscreteVariable(String name, ArrayList<String> values)
	{
		this.name = name;
		pheromone = new ArrayList<Double>();
		this.values = new ArrayList<String>();
		this.values = values;
		for(int i = 0 ; i < values.size() ; i++)
			pheromone.add(DEFAULT_PHEROMONE);
		pheromoneSum = DEFAULT_PHEROMONE*values.size();
		minPheromone = DEFAULT_MIN_PHEROMONE;
		this.min = 0;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return the identifier of the variable
	 */
	public String getName()
	{
		return name;
	}

	/**
	 * @param name The identifier of the variable
	 */
	public void setName(String name)
	{
		this.name = name;
	}

	/**
	 * @return The minimum amount of pheromone a value can have
	 */
	public double getMinPheromone()
	{
		return minPheromone;
	}
	
	/**
	 * @param percentage The percentage of the average pheromone per value to be set as the minimum 
	 * pheromone value. Preferably between 0 and 1.
	 */
	public void setMinPheromone(double percentage)
	{
		percentage = percentage < 0 ? 0 : percentage;
		percentage = percentage > 1 ? 1 : percentage;
		minPheromone = percentage*pheromoneSum/pheromone.size();
	}

	/**
	 * @param minPheromone The minimum absolute amount of pheromone a value can have. Must be 
	 * greater than zero.
	 */
	public void setMinPheromoneAbsolute(double minPheromone) 
	{
		this.minPheromone = minPheromone < 0 ? 0 : minPheromone;
	}

	/**
	 * @param pheromone A list with the pheromone levels for each one of the possible values the 
	 * variable can take. Each level represents a relative measure of the probability of that value 
	 * being chosen when generating new ants. The values must be grater than zero.
	 */
	public void setPheromone(double[] pheromone) 
	{
		pheromoneSum = 0;
		for(int i = 0 ; i < pheromone.length ; i++)
		{
			double value = Math.max(pheromone[i], minPheromone);
			this.pheromone.set(i, value);
			pheromoneSum += value;
		}
	}
	
	/**
	 * @return The number of possible values the variable can take
	 */
	public int getValuesCount()
	{
		return pheromone.size();
	}
	
	/**
	 * Adds a given value to the pheromone level of the value with the provided index.
	 * @param value The index of the value in the pheromone array or the integer value of the 
	 * variable when actual numbers are being used
	 * @param pheromone The amount of pheromone to deposit
	 */
	public void depositPheromone(int value, double pheromone)
	{
		int index = values == null ? value - min : value;
		this.pheromone.set(index, this.pheromone.get(index) + pheromone);
		pheromoneSum += pheromone;
	}
	
	/**
	 * Reduces the levels of pheromone with a uniform reduction for all values. The reduction is 
	 * computed as a percentage of the total pheromone divided into the number of values. No 
	 * pheromone levels are reduced beyond the minimum pheromone value.
	 * @param percentage The fraction of the total pheromone of the variable to be reduced. 
	 * Preferably between 0 and 1.
	 * @return The actual absolute amount of pheromone removed. It may be smaller than the 
	 * requested percentage since no pheromone levels can be reduced past the minimum value.
	 */
	public double evaporatePheromone(double percentage)
	{
		percentage = percentage < 0 ? 0 : percentage;
		percentage = percentage > 1 ? 1 : percentage;
		double total = pheromoneSum*percentage;
		return evaporatePheromoneConstant(total);
	}
	
	/**
	 * Reduces the levels of pheromone with a uniform reduction for all values. The reduction is 
	 * computed from the provided total and dividing it into the number of values. No pheromone 
	 * levels are reduced beyond the minimum pheromone value.
	 * @param total The total amount of pheromone to remove
	 * @return The actual amount of pheromone removed. It may be different from the provided total 
	 * since no pheromone levels can be reduced past the minimum value.
	 */
	public double evaporatePheromoneConstant(double total)
	{
		double ration = total/pheromone.size();
		double sum = 0;
		for(int i = 0 ; i < pheromone.size() ; i++)
		{
			double value = pheromone.get(i);
			double reduction = Math.min(ration, value - minPheromone);
			pheromone.set(i, value - reduction);
			sum += reduction;
		}
		pheromoneSum -= sum;
		return sum;
	}
	
	/**
	 * Generates a value for the variable according to the pheromone level distribution. The 
	 * probability of selecting a given value is proportional to its pheromone level. 
	 * @return The generated value
	 */
	public int generateValue()
	{
		double random = Math.random();
		int index = 0;
		double sum = pheromone.get(index);
		while(sum/pheromoneSum <= random)
		{
			index++;
			sum += pheromone.get(index);
		}
		return values == null ? index + min : index;
	}
	
	/**
	 * Returns the identifier of the specified value of the variable or the value itself if actual 
	 * numbers are being used
	 * @param value The index of the value in the pheromone array or the integer value of the 
	 * variable when actual numbers are being used
	 * @return The identifier of the specified value of the variable or the value itself if actual 
	 * numbers are being used
	 */
	public String getIdentifier(int value)
	{
		return values == null ? String.valueOf(value) : values.get(value);
	}
	
	/**
	 * Resets the pheromone levels to the default values
	 */
	public void reset()
	{
		for(int i = 0 ; i < pheromone.size() ; i++)
			pheromone.set(i, DEFAULT_PHEROMONE);
		pheromoneSum = DEFAULT_PHEROMONE*pheromone.size();
	}
	
}
