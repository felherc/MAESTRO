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
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.TreeSet;

import maestro.ICLogEntry;
import maestro.MAESTROptimizer;
import maestro.gen.GenWrapper;
import maestro.solution.SolutionWrapper;
import utilities.ReportFile;

/**
 * This class writes report files for the execution of MAESTRO
 * @author Felipe Hernández
 */
public class Reports 
{

	// --------------------------------------------------------------------------------------------  
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Report file header 1
	 */
	private final static String REPORT_HEADER_1 = "MAESTRO execution report";
	
	/**
	 * Label for the time the execution of MAESTRO started 
	 */
	private final static String STARTING_TIME_LABEL = "Execution started: ";
	
	/**
	 * Label for the time of the report
	 */
	private final static String TERMINATION_TIME_LABEL = "Execution completed: ";
	
	/**
	 * Label for the number of solutions processed
	 */
	private final static String SOLUTIONS_LABEL = "Solutions processed: ";
	
	/**
	 * Label for the reason for termination
	 */
	private final static String TERMINATION_LABEL = "Terminated: ";

	/**
	 * Identifier of the MAESTRO parameters table
	 */
	private final static String PARAMETERS_TABLE_HEADER = "[MAESTRO parameters]";
	
	/**
	 * Identifier of the <code>updateDist</code> parameter
	 */
	private final static String PARAM_UPDATE_DIST = "Update distribution = ";
	
	/**
	 * Identifier of the <code>genRatio</code> parameter
	 */
	private final static String PARAM_GEN_RATIO = "Generation ratio = ";
	
	/**
	 * Identifier of the <code>genMin</code> parameter
	 */
	private final static String PARAM_GEN_MIN = "Generation minimum = ";
	
	/**
	 * Identifier of the <code>weightPop</code> parameter of the <code>GenManager</code>
	 */
	private final static String PARAM_GEN_WEIGHT_POP = "Population weight (for generation) = ";
	
	/**
	 * Identifier of the <code>weightPart1</code> parameter of the <code>GenManager</code>
	 */
	private final static String PARAM_GEN_WEIGHT_PART1 = "Partition 1 weight (for generation) = ";
	
	/**
	 * Identifier of the <code>weightElite</code> parameter of the <code>GenManager</code>
	 */
	private final static String PARAM_GEN_WEIGHT_ELITE = "Elite set weight (for generation) = ";
	
	/**
	 * Identifier of the <code>popGrowthRate</code> parameter
	 */
	private final static String PARAM_POP_GROWTH_RATE = "Population growth rate = ";
	
	/**
	 * Identifier of the <code>partitionCount</code> parameter
	 */
	private final static String PARAM_PARTITION_COUNT = "Population partitions = ";
	
	/**
	 * Identifier of the <code>minRarityStd</code> parameter
	 */
	private final static String PARAM_MIN_RARITY_STD = "Minimum rarity standard = ";
	
	/**
	 * Identifier of the <code>maxRarityStd</code> parameter
	 */
	private final static String PARAM_MAX_RARITY_STD = "Maximum rarity standard = ";
	
	/**
	 * Identifier of the <code>ictRefuce</code> parameter
	 */
	private final static String PARAM_ICT_REFUSE = "Inner cycle termination refuse = ";
	
	/**
	 * Identifier of the <code>ictUniformity</code> parameter
	 */
	private final static String PARAM_ICT_UNIFORMITY = "Inner cycle termination uniformity = ";
	
	/**
	 * Identifier of the <code>restartElitism</code> parameter
	 */
	private final static String PARAM_RESTART_ELITISM = "Restart elitism method = ";
	
	/**
	 * Identifier of the <code>threadCount</code> parameter
	 */
	private final static String PARAM_THREAD_COUNT = "Number of threads = ";
	
	/**
	 * Identifier of the <code>timeLimit</code> parameter
	 */
	private final static String PARAM_TIME_LIMIT = "Time limit = ";
	
	/**
	 * Identifier of the <code>solutionLimit</code> parameter
	 */
	private final static String PARAM_SOLUTION_LIMIT = "Solution limit = ";
	
	/**
	 * <code>restartElitism</code> <b>no</b> option
	 */
	private final static String RESTART_ELITISM_NO = "no";
	
	/**
	 * <code>restartElitism</code> <b>always</b> option
	 */
	private final static String RESTART_ELITISM_ALWAYS = "always";
	
	/**
	 * <code>restartElitism</code> <b>depends</b> option
	 */
	private final static String RESTART_ELITISM_DEPENDS = "depends";
	
	/**
	 * Identifier of the generators table
	 */
	private final static String GENERATORS_TABLE_HEADER = "[Generator methods]";
	
	/**
	 * Header for the generator method column
	 */
	public final static String GENERATORS_ID = "Generator";
	
	/**
	 * Header for the number of total solutions generated column
	 */
	private final static String GEN_TOTAL_ID = "Total solutions";
	
	/**
	 * Header for the generator parameters column
	 */
	private final static String PARAMETERS_ID = "Parameters";
	
	/**
	 * Identifier of the solution generation history table
	 */
	private final static String GEN_HIST_TABLE_HEADER = "[Generator method use]";
	
	/**
	 * Header for the generation column
	 */
	private final static String GENERATION_ID = "Generation";
	
	/**
	 * Header for the number of solutions generated column
	 */
	private final static String GEN_SOLUTIONS_ID = "Solutions generated";
	
	/**
	 * Header for the time of generation column
	 */
	private final static String GEN_TIME = "Total time (ms)";
	
	/**
	 * Header for the time of generation per solution column
	 */
	private final static String GEN_TIME_PER_SOLUTION = "Time for each (ms)";
	
	/**
	 * Identifier of the inner cycles table
	 */
	private final static String IC_HIST_TABLE_HEADER = "[Inner cycles]";
	
	/**
	 * Header for the inner cycle number column
	 */
	public final static String INNER_CYCLE_ID = "Inner cycle";
	
	/**
	 * Header for the population size column
	 */
	private final static String POPULATION_ID = "Population size";
	
	/**
	 * Header for the number of solutions processed column
	 */
	private final static String SOLUTION_COUNT_ID = "Solutions processed";
	
	/**
	 * Header for the found better solution column
	 */
	private final static String FOUND_BETTER_ID = "Found better solution";
	
	/**
	 * Header for the inner cycle termination column
	 */
	private final static String IC_TERMINATION_ID = "Termination reason";
	
	/**
	 * Identifier of refuse inner cycle termination
	 */
	private final static String IC_TERM_REFUSE = "Refuse";
	
	/**
	 * Identifier of uniformity inner cycle termination
	 */
	private final static String IC_TERM_UNIFORMITY = "Uniformity";
	
	/**
	 * Identifier of outer cycle ended inner cycle termination
	 */
	private final static String IC_TERM_OC_ENDED = "Optimization ended";
	
	/**
	 * Identifier of the best solutions table
	 */
	private final static String BEST_SOLUTIONS_TABLE_HEADER = "[Best solutions]";
	
	/**
	 * Identifier of the elite solution history table
	 */
	private final static String ELITE_SOLUTIONS_TABLE_HEADER = "[Elite solution history]";
	
	/**
	 * Identifier of the table with all the solutions
	 */
	private final static String ALL_SOLUTIONS_TABLE_HEADER = "[All solutions]";
	
	/**
	 * Header for the solutions' names column
	 */
	public final static String SOLUTIONS_ID = "Solution";
	
	/**
	 * Identifier for the discrete variables
	 */
	private final static String DISC_VARIABLES_ID = "Discrete variables";
	
	/**
	 * Identifier for the continuous variables
	 */
	private final static String CONT_VARIABLES_ID = "Continuous variables";
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Writes a report file (or files) containing the results from the execution of MAESTRO. If 
	 * there are too many lines in the resulting file, the report is written to a various files 
	 * with a consecutive number after the provided name.
	 * @param manager The MEASTRO manager that ran the optimization process
	 * @param fileRoute The system route and file name for the report to be written
	 * @param writeConfig True if the parameters and generators of MAESTRO should be included
	 * @param writeICHist True if the history of inner cycles should be included
	 * @param writeGenHist True if the history of generators used should be included
	 * @param writeEliteHist True if the history of best-so-far solutions should be included
	 * @param writeAllSolutions True if all solutions should be included
	 * @return A list with the routes and names of the report files
	 * @throws IOException If the file (or files) can not be created
	 */
	public static ArrayList<String> writeReport(MAESTROptimizer manager, String fileRoute, 
							boolean writeConfig, boolean writeICHist, boolean writeGenHist,
							boolean writeEliteHist, boolean writeAllSolutions) throws IOException
	{
		ReportFile file = new ReportFile();
		
		// Write header
		file.addLineToHeader(REPORT_HEADER_1);
		DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy K:mm:ss a");
		Date date = new Date(manager.getStartTime());
		file.addLineToHeader(STARTING_TIME_LABEL + dateFormat.format(date));
        date = new Date();
        file.addLineToHeader(TERMINATION_TIME_LABEL + dateFormat.format(date));
        file.addLineToHeader(SOLUTIONS_LABEL + manager.getSolutionCount());
        file.addLineToHeader(TERMINATION_LABEL + manager.getTerminationMessage());
        file.addLineToHeader("");
        
        if(writeConfig)
        {
        	String restartElitism =  manager.getRestartElitism() == 
        			MAESTROptimizer.RESTART_ELITISM_NO ? 		RESTART_ELITISM_NO : 
	        		(manager.getRestartElitism() == 
	        		MAESTROptimizer.RESTART_ELITISM_ALWAYS ?	RESTART_ELITISM_ALWAYS :
	        		(manager.getRestartElitism() == 
	        		MAESTROptimizer.RESTART_ELITISM_DEPENDS ?	RESTART_ELITISM_DEPENDS : "N/A"));
        	
        	// Add general parameters
	        int paramTable = file.addTable();
	        file.addLineToTableHeader(paramTable, PARAMETERS_TABLE_HEADER);
	        file.addLineToTableBody(paramTable, PARAM_UPDATE_DIST + "\t" 
	        										+ manager.getUpdateDist());
	        file.addLineToTableBody(paramTable, PARAM_GEN_RATIO + "\t" 
													+ manager.getGenRatio());
	        file.addLineToTableBody(paramTable, PARAM_GEN_MIN + "\t" 
													+ manager.getGenMin());
	        file.addLineToTableBody(paramTable, PARAM_GEN_WEIGHT_POP + "\t" 
													+ manager.getGenWeightPop());
	        file.addLineToTableBody(paramTable, PARAM_GEN_WEIGHT_PART1 + "\t" 
													+ manager.getGenWeightPart1());
	        file.addLineToTableBody(paramTable, PARAM_GEN_WEIGHT_ELITE + "\t" 
													+ manager.getGenWeightElite());
	        file.addLineToTableBody(paramTable, PARAM_POP_GROWTH_RATE + "\t" 
													+ manager.getPopGrowthRate());
	        file.addLineToTableBody(paramTable, PARAM_PARTITION_COUNT + "\t" 
	        										+ manager.getPartitionCount());
	        file.addLineToTableBody(paramTable, PARAM_MIN_RARITY_STD + "\t" 
													+ manager.getMinRarityStd());
	        file.addLineToTableBody(paramTable, PARAM_MAX_RARITY_STD + "\t" 
													+ manager.getMaxRarityStd());
	        file.addLineToTableBody(paramTable, PARAM_ICT_REFUSE + "\t" 
													+ manager.getIctRefuse());
	        file.addLineToTableBody(paramTable, PARAM_ICT_UNIFORMITY + "\t" 
													+ manager.getIctUniformity());
	        file.addLineToTableBody(paramTable, PARAM_RESTART_ELITISM + "\t" 
													+ restartElitism);
	        file.addLineToTableBody(paramTable, PARAM_THREAD_COUNT + "\t" 
													+ manager.getThreadCount());
	        file.addLineToTableBody(paramTable, PARAM_TIME_LIMIT + "\t" 
													+ manager.getTimeLimit() + " ms");
	        file.addLineToTableBody(paramTable, PARAM_SOLUTION_LIMIT + "\t" 
													+ manager.getSolutionLimit());
	        file.addLineToTableFooter(paramTable, "");
	        
	        // Add generator configuration
	        int genTable = file.addTable();
	        file.addLineToTableHeader(genTable, GENERATORS_TABLE_HEADER);
	        file.addLineToTableHeader(genTable, GENERATORS_ID + "\t" + GEN_TOTAL_ID + "\t" 
	        									+ PARAMETERS_ID);
	        ArrayList<GenWrapper> generators = manager.getGenerators();
	        for(GenWrapper generator : generators)
	        {
	        	String line = "";
	        	line += generator.getId() + " (" + generator.getShortId() + ")";
	        	line += "\t" + generator.getGenTotal();
	        	line += "\t" + generator.getParamSummary();
	        	file.addLineToTableBody(genTable, line);
	        }
	        file.addLineToTableFooter(genTable, "");
        }
        
        // Add inner cycle history
        if(writeICHist)
        {
        	int icHistTable = file.addTable();
        	file.addLineToTableHeader(icHistTable, IC_HIST_TABLE_HEADER);
        	file.addLineToTableHeader(icHistTable, INNER_CYCLE_ID + "\t" + POPULATION_ID + "\t" 
        				+ SOLUTION_COUNT_ID + "\t" + FOUND_BETTER_ID + "\t" + IC_TERMINATION_ID);
        	ArrayList<ICLogEntry> icLog = manager.getICLog();
        	for(int i = 0 ; i < icLog.size() ; i++)
        	{
        		ICLogEntry entry = icLog.get(i);
        		String line = (i + 1) + "";
        		line += "\t" + entry.getPopulation();
        		line += "\t" + entry.getSolutions();
        		line += "\t" + entry.foundBest();
        		int term = entry.getTermination();
        		line += "\t" + (term == ICLogEntry.REFUSED ? IC_TERM_REFUSE : 
        				(term == ICLogEntry.UNIFORMITY ? IC_TERM_UNIFORMITY : IC_TERM_OC_ENDED));
        		file.addLineToTableBody(icHistTable, line);
        	}
        	file.addLineToTableFooter(icHistTable, "");
        }
        
    	// Add generator use history
        if(writeGenHist)
        {
        	int genHistTable = file.addTable();
        	file.addLineToTableHeader(genHistTable, GEN_HIST_TABLE_HEADER);
        	file.addLineToTableHeader(genHistTable, GENERATION_ID + "\t" + INNER_CYCLE_ID + "\t" 
        										+ GENERATORS_ID + "\t" + GEN_SOLUTIONS_ID + "\t" 
        										+ GEN_TIME + "\t" + GEN_TIME_PER_SOLUTION);
        	ArrayList<String> gens = manager.getGenHist();
        	for(int i = 0 ; i < gens.size() ; i++)
        		file.addLineToTableBody(genHistTable, (i + 1) + "\t" + gens.get(i));
        	file.addLineToTableFooter(genHistTable, "");
        }
        
        // Add best solutions
        TreeSet<SolutionWrapper> bestSolutions = manager.getBestSolutions();
        if(bestSolutions != null)
	        if(bestSolutions.size() > 0)
	        {
	        	int bestTable = file.addTable();
	        	file.addLineToTableHeader(bestTable, BEST_SOLUTIONS_TABLE_HEADER);
	        	addSolutionHeader(manager, bestSolutions.first(), file, bestTable);
	        	Iterator<SolutionWrapper> solutions = bestSolutions.descendingIterator();
	        	addSolutionList(manager, solutions, file, bestTable);
	        	file.addLineToTableFooter(bestTable, "");
	        }
        
        // Add elite solution history
        if(writeEliteHist)
        {
        	ArrayList<SolutionWrapper> eliteHistory = manager.getEliteHistory();
        	if(eliteHistory.size() > 0)
        	{
        		int eliteTable = file.addTable();
        		file.addLineToTableHeader(eliteTable, ELITE_SOLUTIONS_TABLE_HEADER);
        		addSolutionHeader(manager, eliteHistory.get(0), file, eliteTable);
        		Iterator<SolutionWrapper> solutions = eliteHistory.iterator();
	        	addSolutionList(manager, solutions, file, eliteTable);
	        	file.addLineToTableFooter(eliteTable, "");
        	}
        }
        
        // Add all solutions
        if(writeAllSolutions)
        {
        	ArrayList<SolutionWrapper> allSolutions = manager.getAllSolutions();
        	if(allSolutions.size() > 0)
        	{
        		int allTable = file.addTable();
        		file.addLineToTableHeader(allTable, ALL_SOLUTIONS_TABLE_HEADER);
        		addSolutionHeader(manager, allSolutions.get(0), file, allTable);
        		Iterator<SolutionWrapper> solutions = allSolutions.iterator();
	        	addSolutionList(manager, solutions, file, allTable);
	        	file.addLineToTableFooter(allTable, "");
        	}
        }
        
        // Save to file
        return file.writeFile(fileRoute);
	}
	
	/**
	 * Writes the header of the solution table
	 * @param manager The ACO manager that ran the optimization process
	 * @param solution Any solution of the process
	 * @param file The report file object to write to
	 * @param table The table to write to
	 */
	private static void addSolutionHeader(MAESTROptimizer manager, SolutionWrapper solution, 
											ReportFile file, int table)
	{
		// Retrieve variable count
		ArrayList<Integer> discValues = solution.getSolution().getDiscValues();
		ArrayList<Double> contValues = solution.getSolution().getContValues();
		
		int total = 2;
		int discCount = discValues == null ? 0 : discValues.size();
		int contCount = contValues == null ? 0 : contValues.size();
		
		total += discCount;
		total += contCount;
		
		String[] line1 = new String[total];
		String[] line2 = new String[total];
		String fitnessHeader = solution.getSolution().getReportHeader();
		String space = "";
		for(int i = 0 ; i < fitnessHeader.split("\t").length - 1 ; i++)
			space += "\t";
		
		line1[0] = "\t\t";
		line1[1] = space;
		line2[0] = SOLUTIONS_ID + "\t" + INNER_CYCLE_ID + "\t" + GENERATORS_ID;
		line2[1] = fitnessHeader;
		
		int index = 2;
		int variableIndex = 0;
		
		// Add header for discrete variables
		if(discCount > 0)
		{
			variableIndex = 0;
			line1[index] = DISC_VARIABLES_ID;
			line2[index] = manager.getDiscVarName(variableIndex);
			index++;
			variableIndex++;
			int max = index + discCount - 1;
			for(int i = index ; i < max ; i++)
			{
				line1[i] = "";
				line2[i] = manager.getDiscVarName(variableIndex);
				index++;
				variableIndex++;
			}
		}
		
		// Add header for continuous variables
		if(contCount > 0)
		{
			variableIndex = 0;
			line1[index] = CONT_VARIABLES_ID;
			line2[index] = manager.getContVarName(variableIndex);
			index++;
			variableIndex++;
			int max = index + contCount - 1;
			for(int i = index ; i < max ; i++)
			{
				line1[i] = "";
				line2[i] = manager.getContVarName(variableIndex);
				index++;
				variableIndex++;
			}
		}
		
		// Add header lines to report
		file.addLineToTableHeader(table, line1);
		file.addLineToTableHeader(table, line2);
	}
	
	/**
	 * Writes the list of solutions to the specified table in a report file 
	 * @param manager The MAESTRO manager that ran the optimization process
	 * @param solutions The list of the solutions to write
	 * @param file The report file object to write to
	 * @param table The table to write to
	 */
	private static void addSolutionList(MAESTROptimizer manager, 
								Iterator<SolutionWrapper> solutions, ReportFile file, int table)
	{
		while(solutions.hasNext())
		{
			SolutionWrapper solution = solutions.next();
			
			// Retrieve variable count
			ArrayList<Integer> discValues = solution.getSolution().getDiscValues();
			ArrayList<Double> contValues = solution.getSolution().getContValues();
			
			int total = 4;
			int discCount = discValues == null ? 0 : discValues.size();
			int contCount = contValues == null ? 0 : contValues.size();
			
			total += discCount;
			total += contCount;
			
			String[] line = new String[total];
			line[0] = solution.getId();
			line[1] = solution.getIcIndex() + "";
			line[2] = manager.getGeneratorShortId(solution.getGenIndex());
			line[3] = solution.getSolution().getReport();
			int index = 4;
			
			// Add discrete values
			if(discCount > 0)
				for(int i = 0 ; i < discValues.size() ; i++)
				{
					int value = discValues.get(i);
					line[index] = String.valueOf(manager.getDiscValueID(i, value));
					index++;
				}
			
			// Add continuous values
			if(contCount > 0)
				for(int i = 0 ; i < contValues.size() ; i++)
				{
					line[index] = String.valueOf(contValues.get(i));
					index++;
				}
			
			// Add solution line to table
			file.addLineToTableBody(table, line);
		}
	}

}
