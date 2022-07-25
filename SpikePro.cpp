#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include <climits>
#include <queue>
#include <fstream>
#include "edlib.h"
#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include "CSVparser.hpp"
#include <map>
using namespace std;



int readFastaSequences(const char* path, vector< vector<char> >* seqs);

void printAlignment(const char* query, const char* target,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const EdlibAlignMode modeCode);


int addition (char a, char b);



// For debugging
void printSeq(const vector<char> &seq) {
    for (int i = 0; i < static_cast<int>(seq.size()); i++)
        printf("%d ", seq[i]);
    printf("\n");
}

int main(int argc, char * const argv[]) {

    //----------------------------- PARSE COMMAND LINE ------------------------//
    // If true, there will be no output.
    bool silent = false;
    // Alignment mode.
    char mode[16] = "NW";
    // How many best sequences (those with smallest score) do we want.
    // If 0, then we want them all.
    int numBestSeqs = 0;
    bool findAlignment = false;
    bool findStartLocations = false;
    int option;
    int kArg = -1;
    int numRepeats = 1;
    // If "STD" or "EXT", cigar string will be printed. if "NICE" nice representation
    // of alignment will be printed.
    char alignmentFormat[16] = "NICE";

    bool invalidOption = false;
//    while ((option = getopt(argc, argv, "m:n:k:f:r:spl")) >= 0) {
 //       switch (option) {
 //       case 'm': strcpy(mode, optarg); break;
 //       case 'n': numBestSeqs = atoi(optarg); break;
 //       case 'k': kArg = atoi(optarg); break;
//        case 'f': strcpy(alignmentFormat, optarg); break;
//        case 's': silent = true; break;
//        case 'p': findAlignment = true; break;
//        case 'l': findStartLocations = true; break;
//        case 'r': numRepeats = atoi(optarg); break;
//        default: invalidOption = true;
//        }
//    }

           findAlignment =true;


    if (optind + 2 != argc || invalidOption) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: %s <target.fasta> go\n", argv[0]);
        fprintf(stderr, "\n");
        return 1;
    }
    //-------------------------------------------------------------------------//

    if (strcmp(alignmentFormat, "NICE") && strcmp(alignmentFormat, "CIG_STD") &&
        strcmp(alignmentFormat, "CIG_EXT")) {
        printf("Invalid alignment path format (-f)!\n");
        return 1;
    }

    EdlibAlignMode modeCode;
    if (!strcmp(mode, "SHW"))
        modeCode = EDLIB_MODE_SHW;
    else if (!strcmp(mode, "HW"))
        modeCode = EDLIB_MODE_HW;
    else if (!strcmp(mode, "NW"))
        modeCode = EDLIB_MODE_NW;
    else {
        printf("Invalid mode (-m)!\n");
        return 1;
    }
    printf("Using %s alignment mode.\n", mode);

    EdlibAlignTask alignTask = EDLIB_TASK_DISTANCE;
    if (findStartLocations) alignTask = EDLIB_TASK_LOC;
    if (findAlignment) alignTask = EDLIB_TASK_PATH;


    int readResult;
    // Read queries
    char* queriesFilepath = argv[optind];
    vector< vector<char> >* querySequences = new vector< vector<char> >();
    printf("\nReading query...\n");
    readResult = readFastaSequences(queriesFilepath, querySequences);
    if (readResult) {
        printf("Error: There is no file with name %s\n", queriesFilepath);
        delete querySequences;
        return 1;
    }
//    int numQueries = querySequences->size();
    int numQueries = 1;
    int queriesTotalLength = 0;
    for (int i = 0; i < numQueries; i++) {
        queriesTotalLength += (*querySequences)[i].size();
    }
    printf("\nRead query, %d residues total.\n", queriesTotalLength);
    

    // Read target
    char* targetFilepath = NULL;
    vector< vector<char> >* targetSequences = new vector< vector<char> >();
    printf("\nReading spike protein fasta file...\n");
    readResult = readFastaSequences("P0DTC2.fasta", targetSequences);
    if (readResult) {
        printf("Error: There is no file with name %s\n", targetFilepath);
        delete querySequences;
        delete targetSequences;
        return 1;
    }
    char* target = (*targetSequences)[0].data();
    int targetLength = (*targetSequences)[0].size();
    printf("Read target, %d residues.\n", targetLength);


    // ----------------------------- MAIN CALCULATION ----------------------------- //
    printf("\nComparing queries to target...\n");
    int* scores = new int[numQueries];
    int** endLocations = new int*[numQueries];
    int** startLocations = new int*[numQueries];
    int* numLocations = new int[numQueries];
    priority_queue<int> bestScores; // Contains numBestSeqs best scores
    int k = kArg;
    unsigned char* alignment = NULL; int alignmentLength;
    clock_t start = clock();

    if (!findAlignment || silent) {
        printf("0/%d", numQueries);
        fflush(stdout);
    }
    for (int i = 0; i < numQueries; i++) {
        char* query = (*querySequences)[i].data();
        int queryLength = (*querySequences)[i].size();

        // Calculate score
        EdlibAlignResult result;
        for (int rep = 0; rep < numRepeats; rep++) {  // Redundant repetition, for performance measurements.
            result = edlibAlign(query, queryLength, target, targetLength,
                                edlibNewAlignConfig(k, modeCode, alignTask, NULL, 0));
            if (rep < numRepeats - 1) edlibFreeAlignResult(result);
        }

        scores[i] = result.editDistance;
        endLocations[i] = result.endLocations;
        startLocations[i] = result.startLocations;
        numLocations[i] = result.numLocations;
        alignment = result.alignment;
        alignmentLength = result.alignmentLength;

        // If we want only numBestSeqs best sequences, update best scores 
        // and adjust k to largest score.
        if (numBestSeqs > 0) {
            if (scores[i] >= 0) {
                bestScores.push(scores[i]);
                if (static_cast<int>(bestScores.size()) > numBestSeqs) {
                    bestScores.pop();
                }
                if (static_cast<int>(bestScores.size()) == numBestSeqs) {
                    k = bestScores.top() - 1;
                    if (kArg >= 0 && kArg < k)
                        k = kArg;
                }
            }
        }
        
        if (!findAlignment || silent) {
            printf("\r%d/%d", i + 1, numQueries);
            fflush(stdout);
        } else {
            // Print alignment if it was found, use first position
            if (alignment) {
                printf("\n");
	      	if (!strcmp(alignmentFormat, "NICE")) {
                    printAlignment(query, target, alignment, alignmentLength,
                                   *(endLocations[i]),  modeCode);
                } else {
                    printf("Cigar:\n");
                    EdlibCigarFormat cigarFormat = !strcmp(alignmentFormat, "CIG_STD") ?
                        EDLIB_CIGAR_STANDARD : EDLIB_CIGAR_EXTENDED;
                    char* cigar =edlibAlignmentToCigar(alignment, alignmentLength, cigarFormat);
                    if (cigar) {
                        printf("%s\n", cigar);
                        free(cigar);
                    } else {
                        printf("Error while printing cigar!\n");
                    }
                }
            }
        }

        if (alignment)
            free(alignment);
    }

    if (!silent && !findAlignment) {
        int scoreLimit = -1; // Only scores <= then scoreLimit will be printed (we consider -1 as infinity)
        printf("\n");

        if (bestScores.size() > 0) {
            printf("%d best scores:\n", static_cast<int>(bestScores.size()));
            scoreLimit = bestScores.top();
        } else {
            printf("Scores:\n");
        }

        printf("<query number>: <score>, <num_locations>, "
               "[(<start_location_in_target>, <end_location_in_target>)]\n");
        for (int i = 0; i < numQueries; i++) {
            if (scores[i] > -1 && (scoreLimit == -1 || scores[i] <= scoreLimit)) {
                printf("#%d: %d  %d", i, scores[i], numLocations[i]);
                if (numLocations[i] > 0) {
                    printf("  [");
                    for (int j = 0; j < numLocations[i]; j++) {
                        printf(" (");
                        if (startLocations[i]) {
                            printf("%d", *(startLocations[i] + j));
                        } else {
                            printf("?");
                        }
                        printf(", %d)", *(endLocations[i] + j));
                    }
                    printf(" ]");
                }
                printf("\n");
            }
        }

    }


//     printf("coaooooo %d",addition('D','A'));

//     clock_t finish = clock();
//     double cpuTime = static_cast<double>(finish-start)/CLOCKS_PER_SEC;
    //   printf("\nCpu time of searching: %lf\n", cpuTime);
    // ---------------------------------------------------------------------------- //

    // Free allocated space
    for (int i = 0; i < numQueries; i++) {
        free(endLocations[i]);
        if (startLocations[i]) free(startLocations[i]);
    }
    delete[] endLocations;
    delete[] startLocations;
    delete[] numLocations;
    delete querySequences;
    delete targetSequences;
    delete[] scores;
    
    return 0;
}




/** Reads sequences from fasta file.
 * @param [in] path Path to fasta file containing sequences.
 * @param [out] seqs Sequences will be stored here, each sequence as vector of letters.
 * @return 0 if all ok, positive number otherwise.
 */
int readFastaSequences(const char* path, vector< vector<char> >* seqs) {
    seqs->clear();
    
    FILE* file = fopen(path, "r");
    if (file == 0)
        return 1;

    bool inHeader = false;
    bool inSequence = false;
    const int buffSize = 4096;
    int nonscrivi =0;
    char buffer[buffSize];
    char vaiii[buffSize];
    while (!feof(file)) {
        int read = fread(buffer, sizeof(char), buffSize, file);
        int readheader = fread(vaiii, sizeof(char), buffSize, file);
	for (int i = 0; i < read; ++i) {
		char c = buffer[i];
		if(inHeader==true&&nonscrivi<3){
                        if(i==0){printf(" > ");}
	               printf("%c",buffer[i]);}
		if(c=='|'){nonscrivi=nonscrivi+1;}
            if (inHeader) { // I do nothing if in header
                if (c == '\n')
                    inHeader = false;
            } else {
                if (c == '>') {
                    inHeader = true;
                    inSequence = false;
	     	} else {
                    if (c == '\r' || c == '\n')
                        continue;
                  //         break;  
		  // If starting new sequence, initialize it.
                    if (inSequence == false) {
                        inSequence = true;
                        seqs->push_back(vector<char>());
                    }
                   seqs->back().push_back(c);
                   
	    	}
            }
        }
    }

    fclose(file);
    return 0;
}


void printAlignment(const char* query, const char* target,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const EdlibAlignMode modeCode) {
    int tIdx = -1;
    int qIdx = -1;
    int uprim =0;
    int udop =0;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != EDLIB_EDOP_INSERT)
	     tIdx--;
	}
    }
    int mut[alignmentLength];



    for (int start = 0; start < alignmentLength; start += 50) {
        // target
        printf("T: ");
        int startTIdx = -1;

	for (int j = start; j < start + 50 && j < alignmentLength; j++) {
	     	if (alignment[j] == EDLIB_EDOP_INSERT)
                 {printf("-");
			 mut[j]=-1;}
	    else
	        { printf("%c", target[++tIdx]);
	                     mut[j]=0;}
                
       	    if (j == start)
                startTIdx = tIdx;
          }
          printf(" (%d - %d)\n", max(startTIdx+1, 1), tIdx+1);

      
	  // match / mismatch
          printf("   ");
          for (int j = start; j < start + 50 && j < alignmentLength; j++) {
             printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
 	    if (alignment[j] != EDLIB_EDOP_MATCH){
		          mut[j]++;
	    } 
   	}

       printf("\n");

        // query
        printf("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
	    if (alignment[j] == EDLIB_EDOP_DELETE)
	    {printf("-");
	     mut[j]=0;
	    }
	    else
	    { printf("%c", query[++qIdx]);
		    
         	    if(query[qIdx]=='X')
		        mut[j]=0;
                    if(query[qIdx]=='*')
                        mut[j]=0; 
	    }
		    if (j == start)
                    startQIdx = qIdx;
       }
        printf(" (%d - %d)\n", max(startQIdx, 0), qIdx);
 
  
  
//        printf("M: ");
//        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
 //       printf("%d",mut[j]);
 //       }  
  
        printf("\n");
    }


       // Check number of variants

         int nofvar =0;
         int count  =0;
	 float fitness = 1;
         float fit_NA=1;
         float fit_B=1;
         float fit_S=1;
         float bloomace2=0;
         float bloomplasma=0;
         float bloomoderna=0;

        for (int f = 0; f< alignmentLength;f++){
            if(mut[f]>0){nofvar=nofvar+1;};
	      };


        printf("\nTotal number of amino acid variants: %d\n \n",nofvar); 
  
      
	for (int f = 0; f< alignmentLength;f++){

            if (alignment[f] == EDLIB_EDOP_INSERT)
                 uprim++;

            if (alignment[f] == EDLIB_EDOP_DELETE)
		 udop++;   

               	
	
		if(mut[f]>0){ 
	 	    

		        printf("Variant %c%d%c, occurrence in GISAID = ", target[f-uprim],f+1-uprim,query[f-udop]);
                      
                        int offset = addition(target[f-uprim],query[f-udop])-1; 
		
			csv::Parser file = csv::Parser("PIO_8.csv");
                   
                        string temp;
                       
			float Fitness[ nofvar ];

                        float Fit_NA [ nofvar ];

                        float Fit_B [ nofvar ];

                        float Fit_S [ nofvar ];

                        float Fit_ACE2 [ nofvar ];

                        float Fit_PLASMA [ nofvar ];

                        float Fit_MODERNA [ nofvar ];


      //                 temp = file[f*19+offset]["MutRate"];  // display : 1997
        
                         temp = file[(f-uprim)*19+offset]["MutRate"];

			float int_1 = stof(temp);

                        float pollo = int_1;

                         printf("%.2f%%\n",pollo);

                       temp = file[(f-uprim)*19+offset]["Phi"];  // display : 1997
 
			float int_2 = stof(temp);

			Fitness[count]= fitness*int_2;

                        fitness = Fitness[count];

  // TEST               printf("Test %f\n",int_2);

                        temp = file[(f-uprim)*19+offset]["phi_na"];  // display : 1997

                        float int_3 = stof(temp);

                        Fit_NA[count]= fit_NA*int_3;

                        fit_NA = Fit_NA[count];


                        temp = file[(f-uprim)*19+offset]["phi_b"];  // display : 1997

                        float int_4 = stof(temp);

                        Fit_B[count]= fit_B*int_4;

                        fit_B = Fit_B[count];


                        temp = file[(f-uprim)*19+offset]["phi_s"];  // display : 1997

                        float int_5 = stof(temp);

                        Fit_S[count]= fit_S*int_5;

                        fit_S = Fit_S[count];


                        temp = file[(f-uprim)*19+offset]["BloomACE2"];  // display : 1997

                        float int_6 = stof(temp);

                        Fit_ACE2[count]= bloomace2+int_6;

                        bloomace2 = Fit_ACE2[count];


                        temp = file[(f-uprim)*19+offset]["BloomPLASMA"];  // display : 1997

                        float int_7 = stof(temp);

                        Fit_PLASMA[count]= bloomplasma+int_7;

                        bloomplasma = Fit_PLASMA[count];


                        temp = file[(f-uprim)*19+offset]["BloomMODERNA"];  // display : 1997

                        float int_8 = stof(temp);

                        Fit_MODERNA[count]= bloomoderna+int_8;

                        bloomoderna = Fit_MODERNA[count];


	                count=count+1;
	  
	             }

                   }

         float FIT_TOT=1;

      //  for (int f = 0; f< nofvar;f++){FIT_TOT=FIT_TOT*Fitness[f]};	
     
       printf("\nSpike protein predicted fitness: \u03A6 = %.2f (stab=%.2f,ACE2=%.2f,nAb=%.2f)\n", fitness,fit_S,fit_B,fit_NA);

         printf("\nExperimental value for ACE2 binding: %.2f", bloomace2);

         printf("\nExperimental value for viral escape from human plasma (COVID-19 infected patients): %.2f", bloomplasma);

                printf("\nExperimental value for viral escape from human plasma (MODERNA vaccinated patients): %.2f\n \n", bloomoderna);


     //    for (int y=0;y<1273;y++)
     //	       printf("%d and %c\n", uso, query[y]);

}    



int addition (char a, char b)
{

  std::map<char, int> mp;
	 
	 
    mp['A']=1; 
    mp['C']=2; 
    mp['D']=3; 
    mp['E']=4; 
    mp['F']=5; 
    mp['G']=6;     	 
    mp['H']=7;
    mp['I']=8;
    mp['K']=9;
    mp['L']=10;
    mp['M']=11;
    mp['N']=12;
    mp['P']=13;
    mp['Q']=14;
    mp['R']=15;
    mp['S']=16;
    mp['T']=17;
    mp['V']=18;
    mp['W']=19;
    mp['Y']=20;


      int uno= mp[a];
      int due= mp[b];

   if(uno>due)
	   return due;
   else return due-1;




}




















