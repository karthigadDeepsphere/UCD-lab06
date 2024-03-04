import sequenceAnalysis as SA

class findUnique:
    """
    This class reads a file of fasta sequences from STDIN, finds the unique subsequences that occur in a single tRNA
    that set has no members that occur among any of the other tRNA sets. Each of the sets is minimized,
    such that no member of that tRNA set is found as a substring of any other member of that set. Finally, prints out
    aligned subsequences to STDOUT.
    """
    def __init__(self):
        """
        Takes in multiple sequences of tRNA from a FASTA file and finds power set of each tRNA sequence,
        then appends each power set to a list.
        Attributes:
            attr1 (list): List of power sets for each tRNA.
            attr2 (list): List of unique sets for each tRNA.
            attr3 (dict): Dictionary with key: count (0-22) & value: [tRNA header, tRNA sequence]
        """
        self.powerSetList = []  # powerSetList contains all 22 power sets for each tRNA.
        self.uniqueList = []  # Used to save 22 unique substring sets for each tRNA.  
        self.headerSequenceDictionary = {}  # Dictionary to save the header and sequence of each tRNA as a list.
        fastaFile = SA.FastAreader()  # Instantiate object of FastAreader class so we can read in file.
        count = 0

        for header, sequence in fastaFile.readFasta():  # Reads FastA file and yields the tRNA header/sequence.
           
            filteredSequence = self.removeCharacters(sequence)  # Removes dashes and underscores from sequence.
            self.headerSequenceDictionary[count] = [header, filteredSequence] # Initialize dictionary values with list = [header, squence].
            mypowerSet = self.powerSet(filteredSequence)  # Finds the power set for each tRNA sequence.
            self.powerSetList.append(mypowerSet)  # Add power set for each tRNA to a list.
            count += 1

    def removeCharacters(self, sequence):
        """
        Removes dashes and underscores from tRNA sequences.
        """
        noDashSequence = sequence.replace('-','')
        noUnderscoreSequence = noDashSequence.replace('_','')
        return noUnderscoreSequence
        

    def powerSet(self, sequence):
        """
        This method returns a set of substrings(power set) for each tRNA sequence.
        Ex: 'ABABC' = ['A', 'B', 'C',
                      'AB', 'BC'
                      'ABA', 'BAB', 'ABC',
                      'ABAB', 'BABC',
                      'ABABC']
        """
        powerSets = set()
        for index in range(len(sequence)):
            size = len(sequence)
            while size > index:
                powerSets.add(sequence[index:size])  # Takes big chunk then smaller chunks and moves on to next letter.
                size -= 1
        return powerSets

    def findUniques(self):
        """
        Computes the union of all other tRNA sets, and then removes that set from the current tRNAset to
        find unique subsequences of each tRNA sequence.
        """

        for eachSet in self.powerSetList:  
            
            bigUnion = set()
            copySet = set()
            copyList = self.powerSetList.copy()  # Create copy of original power set list.
            copySet = eachSet.copy()  # Create copy of current power set in power set list.
            copyList.remove(copySet)  # Remove the current power set from our list.
            for powerSet in copyList:
                bigUnion = bigUnion.union(powerSet)  # The union of all the other tRNA sets.
            copySet.difference_update(bigUnion)  # Update the copySet, removing elements found in bigUnion.
            newSet = copySet.copy() 
            for string1 in copySet:
                uniqueSet = copySet.copy() 
                uniqueSet.remove(string1)
                for string2 in uniqueSet:
                    if string1 in string2:  # Checks if string1 is a substring of string2.
                        if len(string1) < len(string2):  # We only want the minimal form.(A, AB)
                            newSet.discard(string2)  # Gets rid of superstrings.
                      
            self.uniqueList.append(newSet)  # Add the truly uniques sets to a list.
            
    def printSequences(self):
        """
        Prints out information to STDOUT.
        """
        for index in range(0,len(self.headerSequenceDictionary)):
            mySequence = self.headerSequenceDictionary[index]
            header = mySequence[0]
            sequence = mySequence[1]
            
            size = len(sequence)
            print(header)
            print(sequence)
            for position in range(0, size):  # Keeps track of position for each sequence.
                for substring in self.uniqueList[index]: # For each unique subsequence in my unique list.
                    substringLength = len(substring)
                    if substring == sequence[position:position + substringLength]: 
                        output = ('.')*position + substring  # 'position' number of '.''s then append substring.
                        print(output)

def main():
    """
    Creates an object of findUnique class, then reads in a FastA file and finds unique subsequences of each tRNA in FastA file.
    Then prints out the header, sequence and aligned unique substrings with '.' to help reader see alignment.
    """
    mytRNA = findUnique()
    mytRNA.findUniques()
    mytRNA.printSequences()
                 
if __name__ == "__main__":
    main()
