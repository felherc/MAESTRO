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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Queue;
import java.util.TreeSet;
import java.util.concurrent.LinkedBlockingQueue;

import maestro.ContVar;
import maestro.DiscVar;
import maestro.ICLogEntry;
import maestro.Monitor;
import maestro.PopulatorThread;
import maestro.Reports;
import maestro.gen.GenManager;
import maestro.gen.GenWrapper;
import maestro.gen.Generator;
import maestro.pop.IndividualUpdatePopulation;
import maestro.pop.GroupMergePopulation;
import maestro.pop.Partition;
import maestro.pop.Population;
import maestro.solution.Solution;
import maestro.solution.SolutionRoot;
import maestro.solution.SolutionThread;
import maestro.solution.SolutionWrapper;
import probDist.ContProbDist;
import utilities.thread.Executor;
import utilities.thread.ExecutorThread;

/**
 * This class allows to solve combinatorial optimization problems using the MAESTRO algorithm. The 
 * problem may have discrete decision variables, continuous decision variables or both. The fitness
 * of candidate solutions is computed using an object provided by the user, which should be able to 
 * compare pairs of solutions given the values of their decision variables. Any specific 
 * combination of population-based solution-generation methods can be used. By default, a Genetic
 * Algorithm (GA), an Ant Colony Optimization (ACO) algorithm, a Covariance Matrix Adaptation
 * Evolution Strategy (CMA-ES), and a Hill-Climbing (HC) algorithm are used.
 * @author Felipe Hernández
 * @version 2
 */
public class MAESTROptimizer implements Executor
{
	
	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Termination message when the provided Solution object is null
	 */
	public final static String TERMINATION_NO_SOLUTION = "The provided solution is null";
	
	/**
	 * Termination message when no variables have been defined
	 */
	public final static String TERMINATION_NO_VARS = "No optimization variables defined";
	
	/**
	 * Termination message when the time limit has been reached
	 */
	public final static String TERMINATION_TIME = "Time limit reached";
	
	/**
	 * Termination message when the number of solutions to be processed has been reached
	 */
	public final static String TERMINATION_SOLUTIONS = "Solutions count limit reached";
	
	/**
	 * Termination message when a solution meets the user defined criteria
	 */
	public final static String TERMINATION_CONVERGED = "Optimization converged; solution meeting"
															+ " convergence criteria: ";
	
	/**
	 * Identifier of {@link #execute(String)} action: notify termination
	 */
	public final static String NOTIFY_TERMINATION = "Notify termination";
	
	/**
	 * Never copy the best solution found so far after a random restart
	 */
	public final static int RESTART_ELITISM_NO = 0;
	
	/**
	 * Always copy the best solution found so far after a random restart
	 */
	public final static int RESTART_ELITISM_ALWAYS = 1;
	
	/**
	 * Only copy the best solution found so far after a random restart if there was no improvement
	 * on the previous inner cycle
	 */
	public final static int RESTART_ELITISM_DEPENDS = 2;
	
	/**
	 * Default value for the percentage of new solutions accepted compared to the size of the 
	 * population before the distributions of the variables are updated
	 */
	public final static double D_UPDATE_DIST = 0.1;
	
	/**
	 * Default value for {@link #populationUpdating}
	 */
	public final static double D_POPULATION_UPDATING = 1.0;
	
	/**
	 * Default value for the initial number of solutions in the population
	 */
	public final static int D_INIT_POP_SIZE = 25;
	
	/**
	 * Default value for the rate at which the population grows after each inner cycle
	 */
	public final static double D_POP_GROWTH_RATE = 1.5;
	
	/**
	 * Default value for {@link #maxPopSize}
	 */
	public final static int D_MAX_POP_SIZE = 2000;
	
	/**
	 * Default value for the number of partitions in the population
	 */
	public final static int D_PARTITION_COUNT = 4;
	
	/**
	 * Default value for the rarity standard of the second partition of the population
	 */
	public final static double D_MIN_RARITY_STD = 0.1;
	
	/**
	 * Default value for the rarity standard of the last partition of the population
	 */
	public final static double D_MAX_RARITY_STD = 0.7;
	
	/**
	 * Default value for the fraction of the number of solutions in the population that can be 
	 * refused before a new inner cycle starts
	 */
	public final static double D_ICT_REFUSE = 5.0;
	
	/**
	 * Default value for the threshold ratio between the current and initial standard deviation of 
	 * the solutions in the first partition of the population
	 */
	public final static double D_ICT_UNIFORMITY = 1E-10;
	
	/**
	 * Default value for the {@link #restartElitism} parameter
	 */
	public final static int D_RESTART_ELITISM = 2;
	
	/**
	 * Default value for the {@link maestro.v2.gen.GenManager#genRatio} attribute
	 */
	public final static double D_GEN_RATIO = 0.25;
	
	/**
	 * Default value for the {@link maestro.v2.gen.GenManager#genMin} attribute
	 */
	public final static double D_GEN_MIN = 0.05;
	
	/**
	 * Default value for the {@link maestro.v2.gen.GenManager#weightPop} attribute
	 */
	public final static double D_GEN_WEIGHT_POP = 0;
	
	/**
	 * Default value for the {@link maestro.v2.gen.GenManager#weightPart1} attribute
	 */
	public final static double D_GEN_WEIGHT_PART_1 = 0.75;
	
	/**
	 * Default value for the {@link maestro.v2.gen.GenManager#weightElite} attribute
	 */
	public final static double D_GEN_WEIGHT_ELITE = 1.0;
	
	/**
	 * Default value for {@link #eliteHistoryFile}
	 */
	public final static String D_ELITE_HISTORY_FILE = "";
	
	/**
	 * Default value for the number of concurrent threads to process new solutions
	 */
	public final static int D_THREAD_COUNT = 1;

	// --------------------------------------------------------------------------------------------
	// Attributes - Structure
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The initial solution object provided by the user from where to create new solutions
	 */
	private Solution solution;
	
	/**
	 * The monitor object provided by the user to communicate with the caller application
	 */
	private Monitor monitor;
	
	/**
	 * The list of the discrete variables of the problem
	 */
	private ArrayList<DiscVar> discVars;
	
	/**
	 * The list of the continuous variables of the problem
	 */
	private ArrayList<ContVar> contVars;
	
	/**
	 * The object that manages the solution generator algorithms
	 */
	private GenManager genManager;
	
	/**
	 * A queue that stores generated values for the solutions threads to create and analyze new
	 * solutions
	 */
	private Queue<SolutionRoot> genBuffer;
	
	/**
	 * The number of solutions in the current population
	 */
	int popSize;
	
	/**
	 * The list of partitions that make up the current population of solutions. Each partition has
	 * a different rarity standard and are sorted in ascending rarity standard order.
	 */
	private Population population;
	
	/**
	 * An ordered set of the best solutions
	 */
	private Partition bestSolutions;
	
	/**
	 * A "hall of fame" list with the best-so-far solutions found by running MAESTRO. Whenever a 
	 * new best is found, it is added at the end of the list.
	 */
	private ArrayList<SolutionWrapper> eliteHistory;
	
	/**
	 * A list that keeps all of the solutions generated
	 */
	private ArrayList<SolutionWrapper> allSolutions;
	
	/**
	 * Thread object that populates the buffer queue
	 */
	private PopulatorThread populator;
	
	/**
	 * A list with the threads that evaluate the new solutions after being generated. Evaluations 
	 * are run within these threads to take advantage of of multi-core CPUs and speeding up run 
	 * times.
	 */
	private ArrayList<SolutionThread> threads;
	
	// --------------------------------------------------------------------------------------------
	// Attributes - Parameters
	// --------------------------------------------------------------------------------------------
	
	/**
	 * 0.0 if the population should be updated every time a new solution is offered. Non-zero 
	 * values indicate the percentage of solutions to be queued as a percentage of the population
	 * capacity before the population is updated.
	 */
	private double populationUpdating;
	
	/**
	 * The initial number of solutions in the population
	 */
	private int initPopSize;
	
	/**
	 * The rate at which the population grows after each inner cycle
	 */
	private double popGrowthRate;
	
	/**
	 * The maximum number of solutions in the population
	 */
	private int maxPopSize;
	
	/**
	 * The number of partitions in the population
	 */
	private int partitionCount;
	
	/**
	 * The rarity standard for the second partition of the population. The rarity standard for each
	 * partition is computed uniformly between the extreme values from the second and last. The 
	 * first partition has no rarity standard.
	 */
	private double minRarityStd;
	
	/**
	 * The rarity standard for the last partition of the population. The rarity standard for each
	 * partition is computed uniformly between the extreme values from the second and last. The 
	 * first partition has no rarity standard.
	 */
	private double maxRarityStd;
	
	/**
	 * The fraction of the number of solutions in the population that can be refused before a new 
	 * inner cycle starts
	 */
	private double ictRefuse;
	
	/**
	 * The threshold ratio between the current and initial standard deviation of the solutions in 
	 * the first partition of the population. <code>Double.NaN</code> if the uniformity criterion 
	 * is not to be used.
	 */
	private double ictUniformity;
	
	/**
	 * The standard deviation of the values of the solutions in the first partition of the 
	 * population at the beginning of the inner cycle when only random solutions have been
	 * generated
	 */
	private double referenceStDev;
	
	/**
	 * Method for implementing elitism after random restart. Defined by the class constants: <ul>
	 * <li> {@link #RESTART_ELITISM_NO}
	 * <li> {@link #RESTART_ELITISM_ALWAYS}
	 * <li> {@link #RESTART_ELITISM_DEPENDS} </ul>
	 */
	private int restartElitism;
	
	/**
	 * The maximum size of the bestSolutions set
	 */
	private int bestSolutionsSize;
	
	/**
	 * The number of concurrent threads to process new solutions
	 */
	private int threadCount;
	
	/**
	 * The maximum run time in milliseconds
	 */
	private long timeLimit;
	
	/**
	 * The maximum number of solutions to be processed
	 */
	private int solutionLimit;
	
	/**
	 * The file route for the file to save {@link #eliteHistory} to. "" if no file should be 
	 * created.
	 */
	private String eliteHistoryFile;
	
	// --------------------------------------------------------------------------------------------
	// Attributes - Control
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The index of the previous solution added to be processed
	 */
	private int solIndex;
	
	/**
	 * The index of the current inner cycle run, starting from 1
	 */
	private int icIndex;
	
	/**
	 * True if solutions are allowed to be processed. False if new solutions should not be 
	 * processed due to inner cycle termination.
	 */
	private boolean icRunning;
	
	/**
	 * The best solution of the previous inner cycle
	 */
	private SolutionWrapper prevICBest;
	
	/**
	 * The list with the log entries that contains the information of the inner cycle runs
	 */
	private ArrayList<ICLogEntry> icLog;
	
	/**
	 * The optimization starting time in milliseconds
	 */
	private long startTime;
	
	/**
	 * The number of solutions processed so far
	 */
	private int solutionCount;
	
	/**
	 * Termination message for the outer cycle; "" if MAESTRO should continue running
	 */
	private String ocTerminate;
	
	/**
	 * True if the judge has been notified of the termination of the optimization process. False 
	 * otherwise. This flag is used to prevent the notification to be issued more than once.
	 */
	private boolean terminationFlag;
	
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
	public MAESTROptimizer(Solution solution, Monitor monitor, boolean keepHist, int bestCount)
	{		
		this.solution				= solution;
		this.monitor				= monitor;
		discVars					= null;
		contVars					= null;
		genManager = new GenManager();
		genManager.setGenRatio(		D_GEN_RATIO);
		genManager.setGenMin(		D_GEN_MIN);
		genManager.setWeightPop(	D_GEN_WEIGHT_POP);
		genManager.setWeightPart1(	D_GEN_WEIGHT_PART_1);
		genManager.setWeightElite(	D_GEN_WEIGHT_ELITE);
		populationUpdating			= D_POPULATION_UPDATING;
		initPopSize					= D_INIT_POP_SIZE;
		popGrowthRate				= D_POP_GROWTH_RATE;
		maxPopSize					= D_MAX_POP_SIZE;
		eliteHistoryFile			= "";
		solIndex 					= 0;
		partitionCount 				= D_PARTITION_COUNT;
		minRarityStd 				= D_MIN_RARITY_STD;
		maxRarityStd 				= D_MAX_RARITY_STD;
		ictRefuse					= D_ICT_REFUSE;
		ictUniformity				= D_ICT_UNIFORMITY;
		referenceStDev 				= Double.NaN;
		restartElitism				= D_RESTART_ELITISM;
		bestSolutionsSize			= Math.max(1, bestCount);
		allSolutions				= keepHist ? new ArrayList<SolutionWrapper>() : null;
		eliteHistory				= new ArrayList<SolutionWrapper>();
		threadCount					= D_THREAD_COUNT;
		genBuffer					= new LinkedBlockingQueue<SolutionRoot>();
		populator					= null;
		if (solution == null)
		{
			ocTerminate				= TERMINATION_NO_SOLUTION;
			terminate(ocTerminate);
		}
	}
	
	// --------------------------------------------------------------------------------------------
	// Getters and setters
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return The initial solution object provided by the user from where to create new solutions
	 */
	public Solution getSolution() 
	{
		return solution;
	}

	/**
	 * @param solution The initial solution object provided by the user from where to create new 
	 * solutions
	 */
	public void setSolution(Solution solution) 
	{
		this.solution = solution;
	}

	/**
	 * @return The monitor object provided by the user to communicate with the caller application
	 */
	public Monitor getMonitor() 
	{
		return monitor;
	}

	/**
	 * @param monitor The monitor object provided by the user to communicate with the caller 
	 * application
	 */
	public void setMonitor(Monitor monitor) 
	{
		this.monitor = monitor;
	}
	
	/**
	 * Returns the name of the discrete variable
	 * @param index The index of the variable
	 * @return The name of the discrete variable
	 */
	public String getDiscVarName(int index)
	{
		return discVars.get(index).getName();
	}
	
	/**
	 * Returns the identifier of the value of a variable
	 * @param varIndex The index of the variable
	 * @param value The index of the value 
	 * @return The identifier of the value of a variable
	 */
	public String getDiscValueID(int varIndex, int value)
	{
		DiscVar variable = discVars.get(varIndex);
		if(variable.getValues() == null)
			return value + "";
		else
			return variable.getValues().get(value);
	}
	
	/**
	 * Returns the name of the continuous variable
	 * @param index The index of the variable
	 * @return The name of the continuous variable
	 */
	public String getContVarName(int index)
	{
		return contVars.get(index).getName();
	}
	
	/**
	 * @return The percentage of new solutions accepted compared to the size of the population 
	 * before the distributions of the variables are updated
	 */
	public double getUpdateDist()
	{
		return population.getUpdateDist();
	}

	/**
	 * @param updateDist The percentage of new solutions accepted compared to the size of the 
	 * population before the distributions of the variables are updated
	 */
	public void setUpdateDist(double updateDist) 
	{
		population.setUpdateDist(updateDist);
	}
	
	/**
	 * @return The list with the algorithms that generate values for new solutions
	 */
	public ArrayList<GenWrapper> getGenerators()
	{
		return genManager.getGenerators();
	}
	
	/**
	 * Adds a new generator method
	 * @param generator The generator method to be added
	 */
	public void addGenerator(Generator generator)
	{
		genManager.addGenerator(generator);
	}
	
	/**
	 * Adds a new Ant Colony Optimization (ACO) generator method with default parameters
	 */
	public void addACOGenerator()
	{
		genManager.addACOGenerator();
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
		genManager.addACOGenerator(minPheromone, elitism, randomP, q, xi);
	}
	
	/**
	 * Adds a new Metropolis - Ant Colony Optimization (MetroACO) generator method with default 
	 * parameters
	 */
	public void addMetroACOGenerator()
	{
		genManager.addMetroACOGenerator();
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
		genManager.addMetroACOGenerator(greed, independentRatio, uniformProb, groupNumber, 
											groupPercent, bandwidthMethod, bandwidthMult);
	}
	
	/**
	 * Adds a new Genetic Algorithm (GA) generator method with default parameters
	 */
	public void addGAGenerator()
	{
		genManager.addGAGenerator();
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
		genManager.addGAGenerator(selectionMethod, q, kPerc, trunc, points, pointUniform, pUniform, 
									unifMethod, unifDistParam, mutationProb, randomMutation, 
									adjacentMutation, boundaryMutation, gaussianMutation);
	}
	
	/**
	 * Adds a new Hill-Climbing (HC) generator method with default parameters
	 */
	public void addHCGenerator()
	{
		genManager.addHCGenerator();
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
		genManager.addHCGenerator(qSolutions, qPairs, truncSolutions, truncPairs, 
								selectionMethodSolutions, selectionMethodPairs, extent, amplitude);
	}
	
	/**
	 * Adds a new Covariance matrix adaptation evolution strategy (CMA-ES) generator method with 
	 * default parameters
	 */
	public void addCMAESGenerator()
	{
		genManager.addCMAESGenerator();
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
		genManager.addCMAESGenerator(q, trunc);
	}
	
	/**
	 * Adds a new Gradient Descent (GD) generator method with default parameters
	 */
	public void addGDGenerator()
	{
		genManager.addGDGenerator();
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
		genManager.addGDGenerator(greed, solVicinities, powerVicinityWeight, stepSize, amplitude, 
									uniformProb);
	}
	
	/**
	 * Returns the identifier of the generator method
	 * @param genIndex The index of the generator method
	 * @return The identifier of the generator method
	 */
	public String getGeneratorId(int genIndex)
	{
		return genManager.getGeneratorId(genIndex);
	}
	
	/**
	 * Returns the short identifier of the generator method
	 * @param genIndex The index of the generator method
	 * @return The short identifier of the generator method
	 */
	public String getGeneratorShortId(int genIndex)
	{
		return genManager.getGeneratorShortId(genIndex);
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
		return genManager.getGenRatio();
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
		genManager.setGenRatio(genRatio);
	}

	/**
	 * @return The minimum number of solutions that should be generated by each generator as a 
	 * percentage of the population size
	 */
	public double getGenMin() 
	{
		return genManager.getGenMin();
	}

	/**
	 * @param genMin The minimum number of solutions that should be generated by each generator as 
	 * a percentage of the population size
	 */
	public void setGenMin(double genMin) 
	{
		genManager.setGenMin(genMin);
	}
	
	/**
	 * @return The absolute minimum number of solutions that should be generated by each generator
	 */
	public int getAbsGenMin()
	{
		return genManager.getAbsGenMin();
	}

	/**
	 * @param absGenMin The absolute minimum number of solutions that should be generated by each 
	 * generator
	 */
	public void setAbsGenMin(int absGenMin) 
	{
		genManager.setAbsGenMin(absGenMin);
	}
	
	/**
	 * @return The weight given to the percentage of solutions in the entire population that were 
	 * created by a generator to determine its overall relative weight. The overall weight for a 
	 * generator is: <br><i>overallWeight = <code>weightPop</code>*%_population + 
	 * <code>weightPart1</code> *%_first_partition + <code>weightElite</code>*%_elite_history</i>.
	 */
	public double getGenWeightPop()
	{
		return genManager.getWeightPop();
	}

	/**
	 * @param weightPop The weight given to the percentage of solutions in the entire population 
	 * that were created by a generator to determine its overall relative weight. The overall 
	 * weight for a generator is: <br><i>overallWeight = <code>weightPop</code>*%_population + 
	 * <code>weightPart1</code> *%_first_partition + <code>weightElite</code>*%_elite_history</i>.
	 */
	public void setGenWeightPop(double weightPop)
	{
		genManager.setWeightPop(weightPop);
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
		return genManager.getWeightPart1();
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
		genManager.setWeightPart1(weightPart1);
	}

	/**
	 * @return The weight given to the percentage of solutions in the history of elite solutions 
	 * that were created by a generator to determine its overall relative weight. The overall 
	 * weight for a generator is: <br><i>overallWeight = <code>weightPop</code>*%_population + 
	 * <code>weightPart1</code> *%_first_partition + <code>weightElite</code>*%_elite_history</i>.
	 */
	public double getGenWeightElite()
	{
		return genManager.getWeightElite();
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
		genManager.setWeightElite(weightElite);
	}

	/**
	 * @return {@link #populationUpdating}
	 */
	public double getPopulationUpdating() 
	{
		return populationUpdating;
	}

	/**
	 * @param populationUpdating {@link #populationUpdating}
	 */
	public void setPopulationUpdating(double populationUpdating) 
	{
		this.populationUpdating = populationUpdating < 0.0 ? 0.0 : populationUpdating;
	}

	/**
	 * @return The initial size of the population
	 */
	public int getInitPopSize()
	{
		return popSize;
	}
	
	/**
	 * @param initPopSize The initial size of the population. Cannot be smaller than 2.
	 */
	public void setInitPopSize(int initPopSize)
	{
		this.initPopSize = initPopSize < 2 ? D_INIT_POP_SIZE : initPopSize;
	}

	/**
	 * @return The log with the history of the number of generated solutions. Each entry on the log 
	 * is a (Inner cycle identifier, generator short identifier, number of solutions) tuple in a 
	 * string separated by tabs.
	 */
	public ArrayList<String> getGenHist()
	{
		return genManager.getGenHist();
	}
	
	/**
	 * @return The rate at which the population grows after each inner cycle
	 */
	public double getPopGrowthRate() 
	{
		return popGrowthRate;
	}

	/**
	 * @param popGrowthRate The rate at which the population grows after each inner cycle
	 */
	public void setPopGrowthRate(double popGrowthRate) 
	{
		this.popGrowthRate = popGrowthRate <= -1 ? -1 + 1E-12 : popGrowthRate;
	}
	
	/**
	 * @return {@link #maxPopSize}
	 */
	public int getMaxPopSize()
	{
		return maxPopSize;
	}
	
	/**
	 * @param maxPopSize {@link #maxPopSize}
	 */
	public void setMaxPopSize(int maxPopSize)
	{
		this.maxPopSize = maxPopSize > initPopSize ? maxPopSize : initPopSize;
	}

	/**
	 * @return The number of partitions in the population
	 */
	public int getPartitionCount() 
	{
		return partitionCount;
	}

	/**
	 * @param partitionCount The number of partitions in the population
	 */
	public void setPartitionCount(int partitionCount) 
	{
		this.partitionCount = partitionCount < 1 ? 1 : partitionCount;
	}
	
	/**
	 * @return The rarity standard for the second partition of the population. The rarity standard 
	 * for each partition is computed uniformly between the extreme values from the second and 
	 * last. The first partition has no rarity standard.
	 */
	public double getMinRarityStd() 
	{
		return minRarityStd;
	}

	/**
	 * @param minRarityStd The rarity standard for the second partition of the population. The 
	 * rarity standard for each partition is computed uniformly between the extreme values from the 
	 * second and last. The first partition has no rarity standard. The minimum rarity standard
	 * should be greater or equal to 0.0.
	 */
	public void setMinRarityStd(double minRarityStd) 
	{
		this.minRarityStd = minRarityStd < 0 ? 0 : minRarityStd;
		this.minRarityStd = minRarityStd >= 1 ? 1 - 1E-16 : minRarityStd;
	}

	/**
	 * @return The rarity standard for the last partition of the population. The rarity standard 
	 * for each partition is computed uniformly between the extreme values from the second and 
	 * last. The first partition has no rarity standard.
	 */
	public double getMaxRarityStd() 
	{
		return maxRarityStd;
	}

	/**
	 * @param maxRarityStd The rarity standard for the last partition of the population. The rarity 
	 * standard for each partition is computed uniformly between the extreme values from the second 
	 * and last. The first partition has no rarity standard. The maximum rarity standard should be
	 * smaller than 1.0.
	 */
	public void setMaxRarityStd(double maxRarityStd) 
	{
		this.maxRarityStd = maxRarityStd < 0 ? 0 : maxRarityStd;
		this.maxRarityStd = maxRarityStd >= 1 ? 1 - 1E-16 : maxRarityStd;
	}
	
	/**
	 * @return The fraction of the number of solutions in the population that can be refused before 
	 * a new inner cycle starts
	 */
	public double getIctRefuse()
	{
		return ictRefuse;
	}
	
	/**
	 * @param ictRefuse The fraction of the number of solutions in the population that can be 
	 * refused before a new inner cycle starts
	 */
	public void setIctRefuse(double ictRefuse)
	{
		if(Double.isNaN(ictRefuse))
			this.ictRefuse = ictRefuse;
		else
			this.ictRefuse = ictRefuse < 0 ? 0 : ictRefuse;
	}

	/**
	 * @return The threshold ratio between the current and initial standard deviation of the 
	 * solutions in the first partition of the population. <code>Double.NaN</code> if the 
	 * uniformity criterion is not to be used.
	 */
	public double getIctUniformity()
	{
		return ictUniformity;
	}
	
	/**
	 * @param ictUniformity The threshold ratio between the current and initial standard deviation 
	 * of the solutions in the first partition of the population. <code>Double.NaN</code> if the 
	 * uniformity criterion is not to be used.
	 */
	public void setIctUniformity(double ictUniformity)
	{
		if(Double.isNaN(ictUniformity))
			this.ictUniformity = ictUniformity;
		else
		{
			this.ictUniformity = ictUniformity < 0 ? 0 : ictUniformity;
			this.ictUniformity = ictUniformity > 1 ? 1 : ictUniformity;
		}
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
		return restartElitism;
	}

	/**
	 * @param restartElitism The method for implementing elitism after random restart. Defined by 
	 * the class constants: <ul>
	 * <li> {@link #RESTART_ELITISM_NO}
	 * <li> {@link #RESTART_ELITISM_ALWAYS}
	 * <li> {@link #RESTART_ELITISM_DEPENDS} </ul>
	 */
	public void setRestartElitism(int restartElitism)
	{
		this.restartElitism = restartElitism;
	}

	/**
	 * @return {@link #eliteHistoryFile}
	 */
	public String getEliteHistoryFile()
	{
		return eliteHistoryFile;
	}

	/**
	 * @param eliteHistoryFile {@link #eliteHistoryFile}
	 * @throws IOException 
	 */
	public void setEliteHistoryFile(String eliteHistoryFile) throws IOException
	{
		this.eliteHistoryFile = eliteHistoryFile;
		PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(
														eliteHistoryFile, true)));
		out.close();
	}

	/**
	 * @return The list with the log entries that contains the information of the inner cycle runs
	 */
	public ArrayList<ICLogEntry> getICLog()
	{
		return icLog;
	}
	
	/**
	 * @return The maximum size of the bestSolutions set
	 */
	public int getBestSolutionsSize() 
	{
		return bestSolutionsSize;
	}

	/**
	 * @param bestSolutionsSize The maximum size of the bestSolutions set
	 */
	public void setBestSolutionsSize(int bestSolutionsSize) 
	{
		this.bestSolutionsSize = bestSolutionsSize < 1 ? 1 : bestSolutionsSize;
	}

	/**
	 * @return The number of concurrent threads to process new solutions
	 */
	public int getThreadCount() 
	{
		return threadCount;
	}

	/**
	 * @param threadCount The number of concurrent threads to process new solutions. Solutions 
	 * should be able to be processed thread-safely if more than one thread is to be used
	 */
	public void setThreadCount(int threadCount) 
	{
		this.threadCount = threadCount < 1 ? 1 : threadCount;
	}

	/**
	 * @return The maximum run time in milliseconds
	 */
	public long getTimeLimit() 
	{
		return timeLimit;
	}
	
	/**
	 * @return The execution start time in milliseconds
	 */
	public long getStartTime()
	{
		return startTime;
	}

	/**
	 * @param timeLimit The maximum run time in milliseconds
	 */
	public void setTimeLimit(long timeLimit) 
	{
		this.timeLimit = timeLimit < 0 ? 0 : timeLimit;
	}

	/**
	 * @return The maximum number of solutions to be processed
	 */
	public int getSolutionLimit() 
	{
		return solutionLimit;
	}

	/**
	 * @param solutionLimit The maximum number of solutions to be processed
	 */
	public void setSolutionLimit(int solutionLimit) 
	{
		this.solutionLimit = solutionLimit < 0 ? 0 : solutionLimit;
	}
	
	/**
	 * @return The number of solutions processed so far
	 */
	public int getSolutionCount()
	{
		return solutionCount;
	}
	
	/**
	 * @return The best solution found
	 */
	public SolutionWrapper getBestSolution()
	{
		return bestSolutions.getSolutions().last();
	}

	/**
	 * @return An ordered set of the best solutions found
	 */
	public TreeSet<SolutionWrapper> getBestSolutions() 
	{
		return bestSolutions.getSolutions();
	}

	/**
	 * @return A "hall of fame" list with the best-so-far solutions found by running MAESTRO. 
	 * Whenever a new best is found, it is added at the end of the list.
	 */
	public ArrayList<SolutionWrapper> getEliteHistory() 
	{
		return eliteHistory;
	}

	/**
	 * @return A list that keeps all of the solutions generated
	 */
	public ArrayList<SolutionWrapper> getAllSolutions() 
	{
		return allSolutions;
	}
	
	/**
	 * @return The message that explains why the execution of MAESTRO ended. Returns an empty 
	 * string ("") if MAESTRO has not started or if the optimization process is still running.  
	 */
	public String getTerminationMessage()
	{
		return ocTerminate;
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
		DiscVar var = new DiscVar(name, min, count, scalar);
		if(discVars == null)
			discVars = new ArrayList<DiscVar>();
		discVars.add(var);
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
		DiscVar var = new DiscVar(name, values, scalar);
		if(discVars == null)
			discVars = new ArrayList<DiscVar>();
		discVars.add(var);
	}
	
	/**
	 * Adds a new continuous optimization variable
	 * @param name Identifier of the variable
	 * @param min The minimum value the variable can take
	 * @param max The maximum value the variable can take
	 */
	public void addContVar(String name, double min, double max)
	{
		ContVar var = new ContVar(name, min, max);
		if(contVars == null)
			contVars = new ArrayList<ContVar>();
		contVars.add(var);
	}
	
	/**
	 * Validates that the values in the provided array are within the valid range of the 
	 * continuous optimization variables. Any values outside the range are replaced by the 
	 * corresponding limit.
	 * @param contValues The list of values to validate
	 */
	public void validateContValues(ArrayList<Double> contValues)
	{
		if (contValues.size() < contVars.size())
			throw new IllegalArgumentException("Argument has only " + contValues.size() + 
										" values; there are " + contVars.size() + " variables");
		for (int v = 0; v < contVars.size(); v++)
			contValues.set(v, contVars.get(v).validate(contValues.get(v)));
	}
	
	/**
	 * Validates that the values in the provided array are within the valid range of the 
	 * discrete optimization variables. Any values outside the range are replaced by the 
	 * corresponding limit.
	 * @param discValues The list of values to validate
	 */
	public void validateDiscValues(ArrayList<Integer> discValues)
	{
		if (discValues.size() < discVars.size())
			throw new IllegalArgumentException("Argument has only " + discValues.size() + 
										" values; there are " + discVars.size() + " variables");
		for (int v = 0; v < discVars.size(); v++)
			discValues.set(v, discVars.get(v).validate(discValues.get(v)));
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
		if(discVars != null || contVars != null)
		{
			genManager.setVariables(discVars, contVars);
			popSize = initPopSize;
			icIndex = 0;
			prevICBest = null;
			icLog = new ArrayList<ICLogEntry>();
			bestSolutions = new Partition(bestSolutionsSize, 0.0);
			allSolutions = allSolutions == null ? null : new ArrayList<SolutionWrapper>();
			eliteHistory = eliteHistory == null ? null : new ArrayList<SolutionWrapper>();
			threads = new ArrayList<SolutionThread>();
			this.timeLimit = timeLimit > 0 ? timeLimit : 0;
			startTime = System.currentTimeMillis();
			this.solutionLimit = solutionLimit > 0 ? solutionLimit : 0;
			solutionCount = 0;
			ocTerminate = "";
			terminationFlag = false;
			prepareEliteHistoryFile();
			initInnerCycle(false);
		}
		else
		{
			ocTerminate = TERMINATION_NO_VARS;
			terminate(ocTerminate);
		}
	}
	
	/**
	 * Writes the header of the elite history file
	 */
	private void prepareEliteHistoryFile() 
	{
		PrintWriter out;
		try 
		{
			if (!eliteHistoryFile.equals(""))
			{
				out = new PrintWriter(new BufferedWriter(new FileWriter(
																eliteHistoryFile, true)));
				String line		= Reports.SOLUTIONS_ID + "\t" + Reports.INNER_CYCLE_ID 
						+ "\t" + Reports.GENERATORS_ID + "\t";
				line			+= solution.getReportHeader() + "\t";
				if (discVars != null)
					for (int i = 0; i < discVars.size(); i++)
						line		+= getDiscVarName(i) + "\t";
				if (contVars != null)
					for (int i = 0; i < contVars.size(); i++)
						line		+= getContVarName(i) + "\t";
				out.println(line);
				out.close();
			}
		} catch (IOException e) 
		{
			e.printStackTrace();
		}
	}

	/**
	 * Initializes a new inner cycle run by initializing the control variables, the population and
	 * the solution threads
	 * @param popGrowth True if the current population size should be modified for the new inner
	 * cycle run
	 */
	private synchronized void initInnerCycle(boolean popGrowth)
	{
		// Initialize control variables
		icIndex++;
		if (icIndex > 1)
		{
			solIndex -= genBuffer.size();
			genBuffer.clear();
		}
		
		// Stop threads
		for (SolutionThread thread : threads)
			thread.end();
		threads.clear();
		
		// Initialize population
		popSize = popGrowth ? (int)((double)popSize*(1 + popGrowthRate)) : popSize;
		popSize = Math.min(maxPopSize, Math.max(popSize, Math.max(partitionCount, 2)));
		population = populationUpdating <= 0.0
						? new IndividualUpdatePopulation(popSize, partitionCount, minRarityStd, 
															maxRarityStd, discVars, contVars)
						: new GroupMergePopulation(popSize, partitionCount, populationUpdating,
												minRarityStd, maxRarityStd, discVars, contVars);
		icLog.add(new ICLogEntry(popSize));
		
		// Add best solution
		boolean improved = false;
		SolutionWrapper bestSoFar = null;
		if(bestSolutions.getSolutions().size() > 0)
			bestSoFar = bestSolutions.getSolutions().last();
		if(prevICBest != null)
		{
			Partition origBestSoFar = bestSoFar.getPartition();
			Partition origPrevICBest = prevICBest.getPartition();
			Partition temp = new Partition(2, 0.0);
			temp.addSolution(bestSoFar);
			temp.addSolution(prevICBest);
			if(restartElitism == RESTART_ELITISM_ALWAYS)
				population.forceAddSolution(bestSoFar);
			if(bestSoFar.compareTo(prevICBest) > 0)
			{
				improved = true;
				if(restartElitism == RESTART_ELITISM_DEPENDS)
					population.forceAddSolution(bestSoFar);
			}
			bestSoFar.setPartition(origBestSoFar);
			prevICBest.setPartition(origPrevICBest);
		}
		prevICBest = bestSoFar;
		if(icIndex == 1)
			icLog.get(0).setFoundBest(true);
		else if(improved)
			icLog.get(icIndex - 2).setFoundBest(true);
		
		// Add randomly generated solutions
		int randSol = popSize - population.getAddedSinceUpdate();
		randSol = Math.max(randSol, 2*threadCount);
		generateRandomRoots(randSol);
		
		// Initialize generators
		genManager.initInnerCycle(population, popSize);
		
		// Start threads
		for (int i = 0; i < threadCount; i++) 
		{
			SolutionThread thread = new SolutionThread(this, solution);
			threads.add(thread);
		}
		icRunning = true;
		for(SolutionThread thread : threads)
			thread.start();
	}

	/**
	 * Includes a number of randomly-generated roots into the generation buffer
	 * @param solutionCount The number of solutions to generate
	 */
	private void generateRandomRoots(int solutionCount)
	{
		ArrayList<ArrayList<Integer>> discValues = null;
		if(discVars != null)
		{
			discValues = new ArrayList<ArrayList<Integer>>();
			for(int j = 0 ; j < discVars.size() ; j++)
				discValues.add(discVars.get(j).generateRandomValues(solutionCount));			
		}
		ArrayList<ArrayList<Double>> contValues = null;
		if(contVars != null)
		{
			contValues = new ArrayList<ArrayList<Double>>();
			for(int j = 0 ; j < contVars.size() ; j++)
				contValues.add(contVars.get(j).generateRandomValues(solutionCount, true));
		}
		for(int i = 0 ; i < solutionCount ; i++)
		{
			ArrayList<Integer> solDiscVals = discValues == null ? null : new ArrayList<Integer>();
			if(solDiscVals != null)
				for(int j = 0 ; j < discVars.size() ; j++)
					solDiscVals.add(discValues.get(j).get(i));
			ArrayList<Double> solContVals = contValues == null ? null : new ArrayList<Double>();
			if(solContVals != null)
				for(int j = 0 ; j < contVars.size() ; j++)
					solContVals.add(contValues.get(j).get(i));
			SolutionRoot solution = new SolutionRoot(solDiscVals, solContVals);
			solIndex++;
			solution.setIndex(solIndex);
			solution.setIcIndex(icIndex);
			genBuffer.offer(solution);
		}
	}
	
	/**
	 * @return The next solution root to be analyzed
	 */
	public synchronized SolutionRoot getSolutionRoot()
	{		
		if(ocTerminate != "")
			return null;
		
		if(icRunning)
		{
			if(genBuffer.size() == 0)
				populateBuffer();
			else if(genBuffer.size() < threads.size())
			{
				boolean okPopulator = populator == null ? true : !populator.isAlive();
				if(okPopulator) 
				{
					populator = new PopulatorThread(this);
					populator.start();
				}
			}
		}		
		return genBuffer.poll();
	}
	
	/**
	 * Adds new solution roots to the buffer to be processed
	 */
	public synchronized void populateBuffer()
	{	
		if(icRunning)
		{
			if(population.size() < GenManager.MIN_POP_SIZE)
			{				
				generateRandomRoots((int)(genManager.getGenRatio()*population.getCapacity()));
				return;
			}
				
			ArrayList<SolutionRoot> roots = genManager.generateSolutions(icIndex, eliteHistory);
			if(roots != null)
				for(SolutionRoot root : roots)
				{
					solIndex++;
					root.setIndex(solIndex);
					root.setIcIndex(icIndex);
					genBuffer.offer(root);
				}
		}
	}
	
	/**
	 * Adds a predefined solution root to be analyzed and offered to the population
	 * @param root The solution root to add
	 */
	public synchronized void addSolutionRoot(SolutionRoot root)
	{
		solIndex++;
		root.setIndex(solIndex);
		root.setIcIndex(icIndex == 0 ? 1 : icIndex);
		genBuffer.offer(root);
	}
	
	/**
	 * Adds a generated and processed solution to the population. The partition distributions and 
	 * the best solution lists are updated if necessary, and the inner and outer cycle termination 
	 * criteria are checked.
	 * @param solution The solution to be added
	 */
	public synchronized void addSolution(SolutionWrapper solution)
	{
		addToLists(solution);
		
		if (ocTerminate.equals("") && icRunning && solution.getIcIndex() == icIndex)
		{
			// Update population
			population.offerSolution(solution);
			
			// Check for met outer cycle termination criteria
			checkOuterCycleTerminate();
			
			// Check for met inner cycle termination criteria
			if(ocTerminate.equals("") && solution.getGenIndex() != -1)
				checkInnerCycleTerminate(solution);
		}
		
		// Verify for convergence
		if (ocTerminate.equals(""))
			if (solution.getSolution().optimizationConverged())
			{
				ocTerminate = TERMINATION_CONVERGED + solution.getId();
				terminate(ocTerminate);
			}
	}

	/**
	 * Adds a solution to the best solution list, to the elite history and to the list with all the
	 * solutions if applicable.
	 * @param solution The solution to add
	 */
	private void addToLists(SolutionWrapper solution) 
	{
		// Update best solution set
		bestSolutions.addSolution(solution);
		
		// Update best solution history
		if(eliteHistory != null)
		{
			if(eliteHistory.size() == 0)
			{
				eliteHistory.add(solution);
				storeHallOfFamer(solution);
			}
			else
			{
				SolutionWrapper best = eliteHistory.get(eliteHistory.size() - 1);
				if(solution.compareTo(best) > 0)
				{
					eliteHistory.add(solution);
					storeHallOfFamer(solution);
				}
			}
		}
		
		// Add to solution list
		if(allSolutions != null)
			allSolutions.add(solution);
	}
	
	/**
	 * Appends the summary of the provided solution to the elite history file
	 * @param solution
	 */
	private void storeHallOfFamer(SolutionWrapper solution)
	{
		try 
		{
			if (!eliteHistoryFile.equals(""))
			{
				PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(
															eliteHistoryFile, true)));
				String line		= solution.getId() + "\t";
				line			+= solution.getIcIndex() + "\t";
				line			+= getGeneratorShortId(solution.getGenIndex()) + "\t";
				line			+= solution.getSolution().getReport() + "\t";
				ArrayList<Integer> discValues	= solution.getSolution().getDiscValues();
				ArrayList<Double> contValues	= solution.getSolution().getContValues();
				if (discValues != null)
					for (int i = 0 ; i < discValues.size() ; i++)
					{
						int value = discValues.get(i);
						line		+= String.valueOf(getDiscValueID(i, value)) + "\t";
					}
				if (contValues != null)
					for (Double value : contValues)
						line		+= value + "\t";
				out.println(line);
				out.close();
			}
		} catch (IOException e) 
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Verifies if the optimization process should end
	 */
	private void checkOuterCycleTerminate() 
	{
		long time = System.currentTimeMillis();
		if(time - startTime >= timeLimit)
			ocTerminate = TERMINATION_TIME;
		solutionCount++;
		if(solutionCount >= solutionLimit)
			ocTerminate = TERMINATION_SOLUTIONS;
		if(!ocTerminate.equals(""))
			terminate(ocTerminate);
	}
	
	/**
	 * Verifies if the current inner cycle should end
	 * @param solution The solution being added
	 */
	private void checkInnerCycleTerminate(SolutionWrapper solution) 
	{
		// Check for met refuse count criterion
		boolean refuseCriterion		= false;
		int refuseCount				= population.getRefuseCount();
		if (!Double.isNaN(ictRefuse))
			refuseCriterion			= refuseCount > (int)(ictRefuse*(double)popSize);
			
		// Check for met uniformity criterion
		boolean uniformityCriterion	= false;
		if (!Double.isNaN(referenceStDev) && refuseCount == 0)
			uniformityCriterion = population.computeFirstPartStDev() <= ictUniformity*referenceStDev;
			 
		if (icRunning)
		{
			if (refuseCriterion || uniformityCriterion)
			{
				// Write current inner cycle log entry and start a new inner cycle
				icRunning			= false;
				int termination		= refuseCriterion ? ICLogEntry.REFUSED 
													: ICLogEntry.UNIFORMITY;
				int accumSolutions	= icLog.size() < 2 ? 0 : 
										icLog.get(icIndex - 2).getAccumSolutions();
				ICLogEntry finished	= icLog.get(icIndex - 1);
				finished.setSolutions(solutionCount	- accumSolutions);
				finished.setAccumSolutions(solutionCount);
				finished.setTermination(termination);
				initInnerCycle(true);
			}
		}
	}
	
	/**
	 * TODO Write comment
	 */
	public void reset()
	{
		// Stop threads
		for(SolutionThread thread : threads)
			thread.end();
		threads.clear();
		
		// TODO Add additional steps and call Monitor 
	}

	/**
	 * Terminates the optimization process if the message is different form an empty string (""). 
	 * In that case, the method <code>terminate()</code> in the monitor is called.
	 * @param message The reason for terminating the algorithm
	 */
	public synchronized void terminate(String message)
	{
		ocTerminate = message;		
		if(message != "")
		{
			// Notify caller in another thread
			ExecutorThread executor = new ExecutorThread(this, NOTIFY_TERMINATION);
			executor.start();
			//notifyTermination();
			
			// Update inner cycle log
			try 
			{
				int accumSolutions = icLog.size() < 2 ? 0 : 
														icLog.get(icIndex - 2).getAccumSolutions();
				boolean improved = prevICBest == null ? true : 
								(bestSolutions.getSolutions().last().compareTo(prevICBest) > 0 ? 
								true : false);
				ICLogEntry finished = icLog.get(icIndex - 1);
				finished.setSolutions(solutionCount - accumSolutions);
				finished.setAccumSolutions(solutionCount);
				finished.setFoundBest(improved);
				finished.setTermination(ICLogEntry.OC_TERMINATED);
			} catch (Exception e) 
			{
				e.printStackTrace();
			}
			
			// Stop threads
			if (threads != null)
			{
				for (SolutionThread thread : threads)
				{
					thread.end();
					thread.interrupt();
				}
				threads = null;
			}
		}
	}
	
	@Override
	public void execute(String processId) 
	{
		if (processId.equals(NOTIFY_TERMINATION))
			notifyTermination();
	}
	
	/**
	 * Notifies the caller that SINOPSYS terminated the execution
	 */
	private synchronized void notifyTermination()
	{
		if(!terminationFlag)
		{
			if(monitor != null)
				monitor.terminate(ocTerminate);
			terminationFlag = true;
		}
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
		boolean writeEliteHist = eliteHistory == null ? false : true;
		boolean writeAllSolutions = allSolutions == null ? false : true;
		return Reports.writeReport(this, fileRoute, writeConfig, writeICHist, writeGenHist,
									writeEliteHist,	writeAllSolutions);
	}
	
}