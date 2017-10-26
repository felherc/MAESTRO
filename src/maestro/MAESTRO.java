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

package maestro;

import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeSet;

import probDist.ContProbDist;
import maestro.ICLogEntry;
import maestro.MAESTROptimizer;
import maestro.Monitor;
import maestro.gen.GenWrapper;
import maestro.gen.Generator;
import maestro.solution.Solution;
import maestro.solution.SolutionRoot;
import maestro.solution.SolutionWrapper;

/**
 * This class allows to solve combinatorial optimization problems using the MAESTRO algorithm. The 
 * problem may have discrete decision variables, continuous decision variables or both. The fitness
 * of candidate solutions is computed using an object provided by the user which should be able to 
 * compare pairs of solutions given the values of their decision variables. Any specific 
 * combination of population-based solution-generation methods can be used. By default, a Genetic
 * Algorithm (GA), an Ant Colony Optimization (ACO) algorithm and a Covariance Matrix Adaptation
 * Evolution Strategy (CMA-ES) are used.
 * @author Felipe Hernández
 * @version 2
 */
public class MAESTRO
{

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The MAESTRO manager object
	 */
	private MAESTROptimizer optimizer;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------

	/**
	 * Instantiates a new MAESTRO optimization object
	 * @param solution An object that implements the <code>Solution</code> interface which helps
	 * creating managing candidate solutions and comparing them to others
	 * @param monitor An object that implements the <code>Monitor</code> interface which allows
	 * the communication with the optimization process
	 * @param keepHist True if all of the solutions generated are to be stored in main memory. 
	 * These can be retrieved after the process has finished its execution with the 
	 * <code>getAllSolutions()</code> method.
	 * @param bestCount The maximum number of solutions to be stored in the best solution list. The 
	 * best solution list can be retrieved after the execution of the algorithm using the 
	 * <code>getBestSolutions()</code> method.
	 */
	public MAESTRO(Solution solution, Monitor monitor, boolean keepHist, int bestCount)
	{
		optimizer = new MAESTROptimizer(solution, monitor, keepHist, bestCount);
	}
	
	// --------------------------------------------------------------------------------------------
	// Getters and setters
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return The initial solution object provided by the user from where to create new solutions
	 */
	public Solution getSolution() 
	{
		return optimizer.getSolution();
	}

	/**
	 * @param solution The initial solution object provided by the user from where to create new 
	 * solutions
	 */
	public void setSolution(Solution solution) 
	{
		optimizer.setSolution(solution);
	}

	/**
	 * @return The monitor object provided by the user to communicate with the caller application
	 */
	public Monitor getMonitor() 
	{
		return optimizer.getMonitor();
	}

	/**
	 * @param monitor The monitor object provided by the user to communicate with the caller 
	 * application
	 */
	public void setMonitor(Monitor monitor) 
	{
		optimizer.setMonitor(monitor);
	}
	
	/**
	 * Returns the name of the discrete variable
	 * @param index The index of the variable
	 * @return The name of the discrete variable
	 */
	public String getDiscVarName(int index)
	{
		return optimizer.getDiscVarName(index);
	}
	
	/**
	 * Returns the identifier of the value of a variable
	 * @param varIndex The index of the variable
	 * @param value The index of the value 
	 * @return The identifier of the value of a variable
	 */
	public String getDiscValueID(int varIndex, int value)
	{
		return optimizer.getDiscValueID(varIndex, value);
	}
	
	/**
	 * Returns the name of the continuous variable
	 * @param index The index of the variable
	 * @return The name of the continuous variable
	 */
	public String getContVarName(int index)
	{
		return optimizer.getContVarName(index);
	}
	
	/**
	 * @return The percentage of new solutions accepted compared to the size of the population 
	 * before the distributions of the variables are updated
	 */
	public double getUpdateDist()
	{
		return optimizer.getUpdateDist();
	}

	/**
	 * @param updateDist The percentage of new solutions accepted compared to the size of the 
	 * population before the distributions of the variables are updated
	 */
	public void setUpdateDist(double updateDist) 
	{
		optimizer.setUpdateDist(updateDist);
	}
	
	/**
	 * @return The list with the algorithms that generate values for new solutions
	 */
	public ArrayList<GenWrapper> getGenerators()
	{
		return optimizer.getGenerators();
	}
	
	/**
	 * Adds a new generator method
	 * @param generator The generator method to be added
	 */
	public void addGenerator(Generator generator)
	{
		optimizer.addGenerator(generator);
	}
	
	/**
	 * Adds a new Ant Colony Optimization (ACO) generator method with default parameters
	 */
	public void addACOGenerator()
	{
		optimizer.addACOGenerator();
	}
	
	/**
	 * Adds a new Ant Colony Optimization (ACO) generator method
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
	public void addACOGenerator(double minPheromone, double elitism, double randomP, double q, 
								double xi)
	{
		optimizer.addACOGenerator(minPheromone, elitism, randomP, q, xi);
	}
	
	/**
	 * Adds a new Metropolis - Ant Colony Optimization (MetroACO) generator method with default 
	 * parameters
	 */
	public void addMetroACOGenerator()
	{
		optimizer.addMetroACOGenerator();
	}
	
	/**
	 * Adds a new Metropolis - Ant Colony Optimization (MetroACO) generator method
	 * @param greed				{@link maestro.v2.gen.aco.MetroACO#greed}
	 * @param independentRatio	{@link maestro.v2.gen.aco.MetroACO#independentRatio}
	 * @param uniformProb		{@link maestro.v2.gen.aco.MetroACO#uniformProb}
	 * @param groupNumber		{@link maestro.v2.gen.aco.MetroACO#groupNumber}
	 * @param groupPercent		{@link maestro.v2.gen.aco.MetroACO#groupPercent}
	 * @param bandwidthMethod	{@link maestro.v2.gen.aco.MetroACO#bandwidthMethod}
	 * @param bandwidthMult		{@link maestro.v2.gen.aco.MetroACO#bandwidthMult}
	 */
	public void addMetroACOGenerator(double greed, ContProbDist independentRatio, 
										double uniformProb, int groupNumber, double groupPercent, 
										int bandwidthMethod, double bandwidthMult)
	{
		optimizer.addMetroACOGenerator(greed, independentRatio, uniformProb, groupNumber, 
											groupPercent, bandwidthMethod, bandwidthMult);
	}
	
	/**
	 * Adds a new Genetic Algorithm (GA) generator method with default parameters
	 */
	public void addGAGenerator()
	{
		optimizer.addGAGenerator();
	}
	
	/**
	 * Adds a new Genetic Algorithm (GA) generator method
	 * Creates a new instance of the Genetic algorithm generator
	 * @param selectionMethod The index of the selection method as defined by the class constants: 
	 * <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 * @param q This variable defines how the solutions in the population are assigned weights. The 
	 * weights are used to compute the likelihood of being selected. <i>q</i> has to be positive. 
	 * A value of <i>0</i> means only the most fitted solution may be selected. A large value of 
	 * <i>q</i> means that every solution has the same probability of being selected with little 
	 * regard for its ranking. Use <code>Double.NaN</code> if an exact uniform distribution should
	 * be used. <i>q</i> is used in all the selection strategies.
	 * @param kPerc The percentage of solutions to be sampled in the Tournament selection strategy
	 * @param trunc The percentage of solutions that are to be truncated from selection. The 
	 * <i>(1 - <code>trunc</code>)*n</i> best solutions are selected.
	 * @param points The number of splitting points in new solutions for the assignment of data 
	 * from either parent 1, parent 2 or uniform crossover. A value of <i>0</i> means that all the 
	 * values in the new solution are computed using the same method.
	 * @param pointUniform The probability of using uniform crossover on a segment of a new 
	 * solution defined by crossover points. <i>0</i> if only point crossover should be used. 
	 * <i>1</i> if only uniform crossover should be used.
	 * @param pUniform The probability of selecting values from parent 1 in uniform crossover. The 
	 * complement is the probability of selecting values from parent 2.
	 * @param unifMethod The index of the uniform crossover method for continuous variables as 
	 * defined by the class constants: <ul>
	 * <li> <code>UNIF_EITHER_OR</code>: Either-or
	 * <li> <code>UNIF_UDIST</code>: Uniform distribution
	 * <li> <code>UNIF_NORMALDIST</code>: Normal distribution </ul>
	 * @param unifDistParam The parameter for probability distribution uniform crossover methods 
	 * for continuous variables: the probability of the generated value to fall outside the range 
	 * between the values of the parents if the uniform distribution method is used; the standard 
	 * deviation of the normal distributions as a percentage of the range between the values of the 
	 * parents if the normal distribution method is used
	 * @param mutationProb The probability of mutating each value in new solutions after crossover
	 * @param randomMutation The weight of the random mutation operator for scalar discrete 
	 * variables. The probability of using the random mutation operator is its weight over the sum 
	 * of the weights of all the operators. 
	 * @param adjacentMutation The weight of the adjacent mutation operator for scalar discrete 
	 * variables. The probability of using the adjacent mutation operator is its weight over the 
	 * sum of the weights of all the operators.
	 * @param boundaryMutation The weight of the boundary mutation operator for scalar discrete 
	 * variables. The probability of using the boundary mutation operator is its weight over the 
	 * sum of the weights of all the operators.
	 * @param gaussianMutation The percentage of a continuous variable range to use as the standard 
	 * deviation of the normal distribution mutation method for continuous variables. 
	 * <code>Double.NaN</code> if the uniform mutation method should be used instead.
	 */
	public void addGAGenerator(int selectionMethod, double q, double kPerc, double trunc, 
			int points, double pointUniform, double pUniform, int unifMethod, double unifDistParam, 
			double mutationProb, double randomMutation,	double adjacentMutation, 
			double boundaryMutation, double gaussianMutation)
	{
		optimizer.addGAGenerator(selectionMethod, q, kPerc, trunc, points, pointUniform, pUniform, 
								unifMethod, unifDistParam, mutationProb, randomMutation, 
								adjacentMutation, boundaryMutation, gaussianMutation);
	}
	
	/**
	 * Adds a new Hill-Climbing (HC) generator method with default parameters
	 */
	public void addHCGenerator()
	{
		optimizer.addHCGenerator();
	}
	
	/**
	 * Adds a new Hill-Climbing (HC) generator method
	 * @param qSolutions This parameter defines how the solutions in the population are assigned 
	 * weights. The weights are used to compute the likelihood of being selected. <i>qSolutions</i> 
	 * has to be positive. A value of <i>0</i> means only the most fitted solution may be selected. 
	 * A large value of <i>qSolutions</i> means that every solution has the same probability of 
	 * being selected with little regard for its ranking. Use <code>Double.NaN</code> if an exact 
	 * uniform distribution should be used.
	 * @param qPairs This parameter defines how the solution pairs are assigned weights. The 
	 * weights are used to compute the likelihood of being selected. <i>qPairs</i> has to be 
	 * positive. A value of <i>0</i> means only the most fitted pair may be selected. A large value 
	 * of <i>qPairs</i> means that every pair has the same probability of being selected with 
	 * little regard for its ranking. Use <code>Double.NaN</code> if an exact uniform distribution 
	 * should be used.
	 * @param truncSolutions The percentage of the number of solutions that are to be selected for 
	 * constructing pairs. <i>(1 - <code>truncSolutions</code>)*n</i> solutions are selected.
	 * @param truncPairs The percentage of the number of pairs that are to be selected for 
	 * constructing new solutions. The <i>(1 - <code>truncPairs</code>)*n</i> best pairs are 
	 * selected.
	 * @param selectionMethodSolutions The index of the selection method for the solutions as 
	 * defined by the class constants: <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 * @param slectionMethodPairs The index of the selection method for the pairs as defined by the 
	 * class constants: <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 * @param extent The extent parameter of the search range for new solutions. New solutions are 
	 * generated using a normal distribution to randomly sample each of its new values. The mean of 
	 * the distribution is created by projecting the gradient of the high rank parent (<i>p1</i>) and 
	 * the low rank parent (<i>p2</i>) for each variable <i>i</i> using the extent as a percentage
	 * of the difference: <i>mean = p1i + extent*(p1i - p2i)</i>.
	 * @param amplitude The amplitude parameter of the search range for new solutions. New 
	 * solutions are generated using a normal distribution to randomly sample each of its new 
	 * values. The amplitude is the standard deviation of the distribution as a percentage of 
	 * <i>p1i - p2i</i>.
	 */
	public void addHCGenerator(double qSolutions, double qPairs, double truncSolutions, 
								double truncPairs, int selectionMethodSolutions, 
								int selectionMethodPairs, double extent, double amplitude)
	{
		optimizer.addHCGenerator(qSolutions, qPairs, truncSolutions, truncPairs, 
								selectionMethodSolutions, selectionMethodPairs, extent, amplitude);
	}
	
	/**
	 * Adds a new Covariance matrix adaptation evolution strategy (CMA-ES) generator method with 
	 * default parameters
	 */
	public void addCMAESGenerator()
	{
		optimizer.addCMAESGenerator();
	}
	
	/**
	 * Adds a new Covariance matrix adaptation evolution strategy (CMA-ES) generator method
	 * @param q This variable defines how the solutions in the population are assigned weights. The 
	 * weights are used to compute the multivariate normal distribution to generate the new 
	 * solutions from. <i>q</i> has to be positive and bigger than 0. A value close to <i>0</i> 
	 * means only the most fitted solutions may contribute to the distribution. A large value of 
	 * <i>q</i> means that every solution has the same contribution with little regard for its 
	 * ranking. Use <code>Double.NaN</code> if an exact uniform distribution should be used.
	 * @param trunc The percentage of solutions that are to be truncated from selection. The 
	 * <i>(1 - <code>trunc</code>)*n</i> best solutions are selected.
	 */
	public void addCMAESGenerator(double q, double trunc)
	{
		optimizer.addCMAESGenerator(q, trunc);
	}
	
	/**
	 * Adds a new Gradient Descent (GD) generator method with default parameters
	 */
	public void addGDGenerator()
	{
		optimizer.addGDGenerator();
	}
	
	/**
	 * Adds a new Gradient Descent (GD) generator method
	 * @param greed					{@link maestro.v2.gen.gd.GD#greed}
	 * @param solVicinities			{@link maestro.v2.gen.gd.GD#solVicinities}
	 * @param powerVicinityWeight	{@link maestro.v2.gen.gd.GD#powerVicinityWeight}
	 * @param stepSize				{@link maestro.v2.gen.gd.GD#stepSize}
	 * @param amplitude				{@link maestro.v2.gen.gd.GD#amplitude}
	 * @param uniformProb			{@link maestro.v2.gen.gd.GD#uniformProb}
	 */
	public void addGDGenerator(double greed, double solVicinities, double powerVicinityWeight, 
						double stepSize, double amplitude, double uniformProb)
	{
		optimizer.addGDGenerator(greed, solVicinities, powerVicinityWeight, stepSize, amplitude, 
									uniformProb);
	}
	
	/**
	 * Returns the identifier of the generator method
	 * @param genIndex The index of the generator method
	 * @return The identifier of the generator method
	 */
	public String getGeneratorId(int genIndex)
	{
		return optimizer.getGeneratorId(genIndex);
	}
	
	/**
	 * Returns the short identifier of the generator method
	 * @param genIndex The index of the generator method
	 * @return The short identifier of the generator method
	 */
	public String getGeneratorShortId(int genIndex)
	{
		return optimizer.getGeneratorShortId(genIndex);
	}
	
	/**
	 * @return The number of solutions for each generator to generate as a ratio of the population 
	 * size when each is showing the same performance. The number of solutions to generate is given 
	 * by: <i>toGenerate = performanceIndex*<code>genRatio</code>*popSize</i>. The performance 
	 * index is computed as: <i>performanceIndex = 
	 * numberOfGenerators*generatorWeight/totalWeight</i>.
	 */
	public double getGenRatio() 
	{
		return optimizer.getGenRatio();
	}

	/**
	 * @param genRatio The number of solutions for each generator to generate as a ratio of the 
	 * population size when each is showing the same performance. The number of solutions to 
	 * generate is given by: <i>toGenerate = performanceIndex*<code>genRatio</code>*popSize</i>. 
	 * The performance index is computed as: <i>performanceIndex = 
	 * numberOfGenerators*generatorWeight/totalWeight</i>.
	 */
	public void setGenRatio(double genRatio) 
	{
		optimizer.setGenRatio(genRatio);
	}

	/**
	 * @return The minimum percentage of solutions that should be generated by each generator
	 */
	public double getGenMin() 
	{
		return optimizer.getGenMin();
	}

	/**
	 * @param genMin The minimum percentage of solutions that should be generated by each generator
	 */
	public void setGenMin(double genMin) 
	{
		optimizer.setGenMin(genMin);
	}
	
	/**
	 * @return The absolute minimum number of solutions that should be generated by each generator
	 */
	public int getAbsGenMin()
	{
		return optimizer.getAbsGenMin();
	}

	/**
	 * @param absGenMin The absolute minimum number of solutions that should be generated by each 
	 * generator
	 */
	public void setAbsGenMin(int absGenMin) 
	{
		optimizer.setAbsGenMin(absGenMin);
	}
	
	/**
	 * @return The weight given to the percentage of solutions in the entire population that were 
	 * created by a generator to determine its overall relative weight. The overall weight for a 
	 * generator is: <br><i>overallWeight = <code>weightPop</code>*%_population + 
	 * <code>weightPart1</code> *%_first_partition + <code>weightElite</code>*%_elite_history</i>.
	 */
	public double getGenWeightPop()
	{
		return optimizer.getGenWeightPop();
	}

	/**
	 * @param weightPop The weight given to the percentage of solutions in the entire population 
	 * that were created by a generator to determine its overall relative weight. The overall 
	 * weight for a generator is: <br><i>overallWeight = <code>weightPop</code>*%_population + 
	 * <code>weightPart1</code> *%_first_partition + <code>weightElite</code>*%_elite_history</i>.
	 */
	public void setGenWeightPop(double weightPop)
	{
		optimizer.setGenWeightPop(weightPop);
	}

	/**
	 * @return The weight given to the percentage of solutions in the first partition of the 
	 * population that were created by a generator to determine its overall relative weight. The 
	 * overall weight for a generator is: <br><i>overallWeight = <code>weightPop</code>*
	 * %_population + <code>weightPart1</code> *%_first_partition + <code>weightElite</code>*
	 * %_elite_history</i>.
	 */
	public double getGenWeightPart1() 
	{
		return optimizer.getGenWeightPart1();
	}

	/**
	 * @param weightPart1 The weight given to the percentage of solutions in the first partition of 
	 * the population that were created by a generator to determine its overall relative weight. 
	 * The overall weight for a generator is: <br><i>overallWeight = <code>weightPop</code>*
	 * %_population + <code>weightPart1</code> *%_first_partition + <code>weightElite</code>*
	 * %_elite_history</i>.
	 */
	public void setGenWeightPart1(double weightPart1) 
	{
		optimizer.setGenWeightPart1(weightPart1);
	}

	/**
	 * @return The weight given to the percentage of solutions in the history of elite solutions 
	 * that were created by a generator to determine its overall relative weight. The overall 
	 * weight for a generator is: <br><i>overallWeight = <code>weightPop</code>*%_population + 
	 * <code>weightPart1</code> *%_first_partition + <code>weightElite</code>*%_elite_history</i>.
	 */
	public double getGenWeightElite()
	{
		return optimizer.getGenWeightElite();
	}

	/**
	 * @param weightElite The weight given to the percentage of solutions in the history of elite 
	 * solutions that were created by a generator to determine its overall relative weight. The 
	 * overall weight for a generator is: <br><i>overallWeight = <code>weightPop</code>*
	 * %_population + <code>weightPart1</code> *%_first_partition + <code>weightElite</code>*
	 * %_elite_history</i>.
	 */
	public void setGenWeightElite(double weightElite) 
	{
		optimizer.setGenWeightElite(weightElite);
	}
	
	/**
	 * @return 0.0 if the population should be updated every time a new solution is offered. 
	 * Non-zero values indicate the percentage of solutions to be queued as a percentage of the 
	 * population capacity before the population is updated.
	 */
	public double getPopulationUpdating() 
	{
		return optimizer.getPopulationUpdating();
	}

	/**
	 * @param populationUpdating 0.0 if the population should be updated every time a new solution 
	 * is offered. Non-zero values indicate the percentage of solutions to be queued as a 
	 * percentage of the population capacity before the population is updated.
	 */
	public void setPopulationUpdating(double populationUpdating) 
	{
		optimizer.setPopulationUpdating(populationUpdating);
	}
	
	/**
	 * @return The initial size of the population
	 */
	public int getInitPopSize()
	{
		return optimizer.getInitPopSize();
	}
	
	/**
	 * @param initPopSize The initial size of the population. Cannot be smaller than 2.
	 */
	public void setInitPopSize(int initPopSize)
	{
		optimizer.setInitPopSize(initPopSize);
	}

	/**
	 * @return The log with the history of the number of generated solutions. Each entry on the log 
	 * is a (Inner cycle identifier, generator short identifier, number of solutions) tuple in a 
	 * string separated by tabs.
	 */
	public ArrayList<String> getGenHist()
	{
		return optimizer.getGenHist();
	}
	
	/**
	 * @return The rate at which the population growths after each inner cycle
	 */
	public double getPopGrowthRate() 
	{
		return optimizer.getPopGrowthRate();
	}

	/**
	 * @param popGrowthRate The rate at which the population growths after each inner cycle
	 */
	public void setPopGrowthRate(double popGrowthRate) 
	{
		optimizer.setPopGrowthRate(popGrowthRate);
	}
	
	/**
	 * @return The maximum number of solutions in the population
	 */
	public int getMaxPopSize()
	{
		return optimizer.getMaxPopSize();
	}
	
	/**
	 * @param maxPopSize The maximum number of solutions in the population
	 */
	public void setMaxPopSize(int maxPopSize)
	{
		optimizer.setMaxPopSize(maxPopSize);
	}

	/**
	 * @return The number of partitions in the population
	 */
	public int getPartitionCount() 
	{
		return optimizer.getPartitionCount();
	}

	/**
	 * @param partitionCount The number of partitions in the population
	 */
	public void setPartitionCount(int partitionCount) 
	{
		optimizer.setPartitionCount(partitionCount);
	}
	
	/**
	 * @return The rarity standard for the second partition of the population. The rarity standard 
	 * for each partition is computed uniformly between the extreme values from the second and 
	 * last. The first partition has no rarity standard.
	 */
	public double getMinRarityStd() 
	{
		return optimizer.getMinRarityStd();
	}

	/**
	 * @param minRarityStd The rarity standard for the second partition of the population. The 
	 * rarity standard for each partition is computed uniformly between the extreme values from the 
	 * second and last. The first partition has no rarity standard. The minimum rarity standard
	 * should be greater or equal to 0.0.
	 */
	public void setMinRarityStd(double minRarityStd) 
	{
		optimizer.setMinRarityStd(minRarityStd);
	}

	/**
	 * @return The rarity standard for the last partition of the population. The rarity standard 
	 * for each partition is computed uniformly between the extreme values from the second and 
	 * last. The first partition has no rarity standard.
	 */
	public double getMaxRarityStd() 
	{
		return optimizer.getMaxRarityStd();
	}

	/**
	 * @param maxRarityStd The rarity standard for the last partition of the population. The rarity 
	 * standard for each partition is computed uniformly between the extreme values from the second 
	 * and last. The first partition has no rarity standard. The maximum rarity standard should be
	 * smaller than 1.0.
	 */
	public void setMaxRarityStd(double maxRarityStd) 
	{
		optimizer.setMaxRarityStd(maxRarityStd);
	}
	
	/**
	 * @return The fraction of the number of solutions in the population that can be refused before 
	 * a new inner cycle starts
	 */
	public double getIctRefuse()
	{
		return optimizer.getIctRefuse();
	}
	
	/**
	 * @param ictRefuse The fraction of the number of solutions in the population that can be 
	 * refused before a new inner cycle starts
	 */
	public void setIctRefuse(double ictRefuse)
	{
		optimizer.setIctRefuse(ictRefuse);
	}

	/**
	 * @return The threshold ratio between the current and initial standard deviation of the 
	 * solutions in the first partition of the population. <code>Double.NaN</code> if the 
	 * uniformity criterion is not to be used.
	 */
	public double getIctUniformity()
	{
		return optimizer.getIctUniformity();
	}
	
	/**
	 * @param ictUniformity The threshold ratio between the current and initial standard deviation 
	 * of the solutions in the first partition of the population. <code>Double.NaN</code> if the 
	 * uniformity criterion is not to be used.
	 */
	public void setIctUniformity(double ictUniformity)
	{
		optimizer.setIctUniformity(ictUniformity);
	}
	
	/**
	 * @return The method for implementing elitism after random restart. Defined by the class 
	 * constants: <ul>
	 * <li> {@link #RESTART_ELITISM_NO}
	 * <li> {@link #RESTART_ELITISM_ALWAYS}
	 * <li> {@link #RESTART_ELITISM_DEPENDS} </ul>
	 */
	public int getRestartElitism() 
	{
		return optimizer.getRestartElitism();
	}

	/**
	 * @param restartElitism The method for implementing elitism after random restart. Defined by 
	 * the class constants: <ul>
	 * <li> {@link MAESTROptimizer#RESTART_ELITISM_NO}
	 * <li> {@link MAESTROptimizer#RESTART_ELITISM_ALWAYS}
	 * <li> {@link MAESTROptimizer#RESTART_ELITISM_DEPENDS} </ul>
	 */
	public void setRestartElitism(int restartElitism)
	{
		optimizer.setRestartElitism(restartElitism);
	}
	
	/**
	 * @return The file route for the file to save the "hall of fame" list with the best-so-far 
	 * solutions found by running MAESTRO. Whenever a new best is found, it is added at the end of 
	 * the list. to. "" if no file should be created.
	 */
	public String getEliteHistoryFile()
	{
		return optimizer.getEliteHistoryFile();
	}

	/**
	 * @param eliteHistoryFile The file route for the file to save the "hall of fame" list with the 
	 * best-so-far solutions found by running MAESTRO. Whenever a new best is found, it is added at 
	 * the end of the list. to. "" if no file should be created.
	 * @throws IOException 
	 */
	public void setEliteHistoryFile(String eliteHistoryFile) throws IOException
	{
		optimizer.setEliteHistoryFile(eliteHistoryFile);
	}
	
	/**
	 * @return The list with the log entries that contains the information of the inner cycle runs
	 */
	public ArrayList<ICLogEntry> getICLog()
	{
		return optimizer.getICLog();
	}
	
	/**
	 * @return The maximum size of the bestSolutions set
	 */
	public int getBestSolutionsSize() 
	{
		return optimizer.getBestSolutionsSize();
	}

	/**
	 * @param bestSolutionsSize The maximum size of the bestSolutions set
	 */
	public void setBestSolutionsSize(int bestSolutionsSize) 
	{
		optimizer.setBestSolutionsSize(bestSolutionsSize);
	}

	/**
	 * @return The number of concurrent threads to process new solutions
	 */
	public int getThreadCount() 
	{
		return optimizer.getThreadCount();
	}

	/**
	 * @param threadCount The number of concurrent threads to process new solutions. Solutions 
	 * should be able to be processed thread-safely if more than one thread is to be used
	 */
	public void setThreadCount(int threadCount) 
	{
		optimizer.setThreadCount(threadCount);
	}

	/**
	 * @return The maximum run time in milliseconds
	 */
	public long getTimeLimit() 
	{
		return optimizer.getTimeLimit();
	}
	
	/**
	 * @param timeLimit The maximum run time in milliseconds
	 */
	public void setTimeLimit(long timeLimit) 
	{
		optimizer.setTimeLimit(timeLimit);
	}

	/**
	 * @return The maximum number of solutions to be processed
	 */
	public int getSolutionLimit() 
	{
		return optimizer.getSolutionLimit();
	}

	/**
	 * @param solutionLimit The maximum number of solutions to be processed
	 */
	public void setSolutionLimit(int solutionLimit) 
	{
		optimizer.setSolutionLimit(solutionLimit);
	}
	
	/**
	 * @return The number of solutions processed so far
	 */
	public int getSolutionCount()
	{
		return optimizer.getSolutionCount();
	}
	
	/**
	 * @return The best solution found
	 */
	public SolutionWrapper getBestSolution()
	{
		return optimizer.getBestSolution();
	}

	/**
	 * @return An ordered set of the best solutions found
	 */
	public TreeSet<SolutionWrapper> getBestSolutions() 
	{
		return optimizer.getBestSolutions();
	}

	/**
	 * @return A "hall of fame" list with the best-so-far solutions found by running MAESTRO. 
	 * Whenever a new best is found, it is added at the end of the list.
	 */
	public ArrayList<SolutionWrapper> getEliteHistory() 
	{
		return optimizer.getEliteHistory();
	}

	/**
	 * @return A list that keeps all of the solutions generated
	 */
	public ArrayList<SolutionWrapper> getAllSolutions() 
	{
		return optimizer.getAllSolutions();
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Adds a new discrete optimization variable with its possible values. The values are the 
	 * specified number of integers starting from the <code>min</code> value and sequentially on.
	 * @param name The identifier of the variable
	 * @param min The minimum integer value the variable can take
	 * @param count The number of values the variable can take
	 * @param scalar True if the values of the variable correspond to a scale
	 */
	public void addDiscVar(String name, int min, int count, boolean scalar)
	{
		optimizer.addDiscVar(name, min, count, scalar);
	}
	
	/**
	 * Adds a new discrete optimization variable with its possible values. The values are integers 
	 * from 0 and on, that map to the set of identifiers.
	 * @param name The identifier of the variable
	 * @param values The set of identifiers of the variable's possible values
	 * @param scalar True if the values of the variable correspond to a scale
	 */
	public void addDiscVar(String name, ArrayList<String> values, boolean scalar)
	{
		optimizer.addDiscVar(name, values, scalar);
	}
	
	/**
	 * Adds a new continuous optimization variable
	 * @param name Identifier of the variable
	 * @param min The minimum value the variable can take
	 * @param max The maximum value the variable can take
	 */
	public void addContVar(String name, double min, double max)
	{
		optimizer.addContVar(name, min, max);
	}
	
	/**
	 * Adds a predefined solution root to be analyzed and offered to the population. Do not call
	 * this method before starting the optimization as the added solutions will not be processed.
	 * @param root The solution root to add
	 */
	public synchronized void addSolutionRoot(SolutionRoot root)
	{
		optimizer.addSolutionRoot(root);
	}
	
	/**
	 * Validates that the values in the provided array are within the valid range of the 
	 * continuous optimization variables. Any values outside the range are replaced by the 
	 * corresponding limit.
	 * @param contValues The list of values to validate
	 */
	public void validateContValues(ArrayList<Double> contValues)
	{
		optimizer.validateContValues(contValues);
	}
	
	/**
	 * Validates that the values in the provided array are within the valid range of the 
	 * discrete optimization variables. Any values outside the range are replaced by the 
	 * corresponding limit.
	 * @param discValues The list of values to validate
	 */
	public void validateDiscValues(ArrayList<Integer> discValues)
	{
		optimizer.validateDiscValues(discValues);
	}
	
	/**
	 * Initializes the optimization process which runs until the time limit is reached, the maximum
	 * number of solutions is generated or the method <code>terminate()</code> is called. When the 
	 * process is finished, the method <code>terminate()</code> in the <code>Monitor</code> is 
	 * called.
	 * @param timeLimit The time limit for the optimization process in milliseconds
	 * @param solutionLimit The maximum number of candidate solutions that should be processed
	 */
	public synchronized void startOptimization(long timeLimit, int solutionLimit)
	{
		optimizer.startOptimization(timeLimit, solutionLimit);
	}
	
	/**
	 * Terminates the optimization process if the message is different form an empty string (""). 
	 * In that case, the method <code>terminate()</code> in the monitor is called.
	 * @param message The reason for terminating the algorithm
	 */
	public synchronized void terminate(String message)
	{
		optimizer.terminate(message);
	}
	
	/**
	 * TODO Write comment
	 */
	public void reset()
	{
		optimizer.reset();
	}
	
	/**
	 * Writes a report file (or files) containing the results from the execution of MAESTRO. If 
	 * there are too many lines in the resulting file, the report is written to a various files 
	 * with a consecutive number after the provided name.
	 * @param fileRoute The system route and file name for the report to be written
	 * @param writeConfig True if the parameters and generators of MAESTRO should be included
	 * @param writeEliteHist True if the history of best-so-far solutions should be included
	 * @param writeICHist True if the history of inner cycles should be included
	 * @param writeGenHist True if the history of generators used should be included
	 * @return A list with the routes and names of the report files
	 * @throws IOException If the file (or files) can not be created
	 */
	public ArrayList<String> writeReport(String fileRoute, boolean writeConfig, 
									boolean writeICHist, boolean writeGenHist) throws IOException
	{
		return optimizer.writeReport(fileRoute, writeConfig, writeICHist, writeGenHist);
	}
	
}
