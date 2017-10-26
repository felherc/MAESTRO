package maestro.gen.aco;

/**
 * This class stores the names, minimum and maximum values of continuous variables for an Ant 
 * colony Optimization problem
 * @author Felipe Hernández
 */
public class ContinuousVariable 
{
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Identifier of the variable
	 */
	private String name;
	
	/**
	 * The minimum value the variable can take
	 */
	private double minValue;
	
	/**
	 * The maximum value the variable can take
	 */
	private double maxValue;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Initializes a new continuous variable
	 * @param name The identifier of the variable
	 * @param minValue The minimum value the variable can take
	 * @param maxValue The maximum value the variable can take
	 * @throws Exception If the maximum value is not greater than the minimum value
	 */
	public ContinuousVariable(String name, double minValue, double maxValue)
	{
		this.name = name;
		this.minValue = minValue < maxValue ? minValue : maxValue;
		this.maxValue = maxValue > minValue ? maxValue : minValue;
	}

	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return The identifier of the variable
	 */
	public String getName()
	{
		return name;
	}

	/**
	 * @return The minimum value the variable can take
	 */
	public double getMinValue()
	{
		return minValue;
	}

	/**
	 * @return The maximum value the variable can take
	 */
	public double getMaxValue()
	{
		return maxValue;
	}
	
}
