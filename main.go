package main

import (
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// main parses the command-line flags and either builds the index or performs a search.
func main() {
	indexMode := flag.Bool("m", false, "Index mode: build suffix array index")
	searchQueryStr := flag.String("s", "", "Search mode: search for sequence")
	fileName := flag.String("f", "genoma.txt", "Genome file name")
	flag.Parse()

	// Read genome file
	data, err := os.ReadFile(*fileName)
	if err != nil {
		fmt.Println("Error reading genome file:", err)
		os.Exit(1)
	}
	genome := strings.TrimSpace(string(data))

	if *indexMode {
		fmt.Println("Building suffix array index using SAIS algorithm...")
		// Encode the genome: add a sentinel 0 and shift characters so that 0 is reserved.
		encoded, alphabetSize := encodeString(genome)
		sa := SAISEntryPoint(encoded, alphabetSize)
		// Save the suffix array to a file (one index per line)
		err = saveIndex("sa.idx", sa)
		if err != nil {
			fmt.Println("Error saving index:", err)
			os.Exit(1)
		}
		fmt.Println("Index built and saved to sa.idx")
	} else if *searchQueryStr != "" {
		// Load the suffix array index
		sa, err := loadIndex("sa.idx")
		if err != nil {
			fmt.Println("Error loading index:", err)
			os.Exit(1)
		}
		fmt.Printf("Searching for sequence: %s\n", *searchQueryStr)
		positions := searchSequence(genome, sa, *searchQueryStr)
		if len(positions) == 0 {
			fmt.Println("Sequence not found.")
		} else {
			fmt.Println("Sequence found at positions:")
			for _, pos := range positions {
				fmt.Printf("%d ", pos)
			}
			fmt.Println()
		}
	} else {
		fmt.Println("Please provide -m to build index or -s <sequence> to search.")
	}
}

// saveIndex writes the suffix array to a file.
func saveIndex(filename string, sa []int) error {
	var lines []string
	for _, num := range sa {
		lines = append(lines, strconv.Itoa(num))
	}
	content := strings.Join(lines, "\n")
	return os.WriteFile(filename, []byte(content), 0644)
}

// loadIndex reads the suffix array from a file.
func loadIndex(filename string) ([]int, error) {
	data, err := os.ReadFile(filename)
	if err != nil {
		return nil, err
	}
	lines := strings.Split(strings.TrimSpace(string(data)), "\n")
	sa := make([]int, len(lines))
	for i, line := range lines {
		num, err := strconv.Atoi(line)
		if err != nil {
			return nil, err
		}
		sa[i] = num
	}
	return sa, nil
}

// encodeString converts the input string into an int slice and appends a sentinel (0).
// Each character is shifted by +1 so that 0 is reserved as the sentinel.
func encodeString(s string) ([]int, int) {
	n := len(s)
	encoded := make([]int, n+1)
	maxVal := 0
	for i, ch := range s {
		encoded[i] = int(ch) + 1
		if encoded[i] > maxVal {
			maxVal = encoded[i]
		}
	}
	encoded[n] = 0
	// The alphabet size is maxVal+1 (since we shifted all characters by 1)
	return encoded, maxVal + 1
}

func SAISEntryPoint(s []int, K int) []int {
	n := len(s)
	return SAIS(s, K, n, make([]int, n), make([]int, n))
}

// SAIS constructs the suffix array for s using the SAIS algorithm.
// s is expected to have a trailing sentinel (0).
func SAIS(s []int, K int, n int, SA []int, lmsNames []int) []int {
	SA = SA[:n]
	for i := range SA {
		SA[i] = -1
	}
	if n == 0 {
		return SA
	}
	if n == 1 {
		SA[0] = 0
		return SA
	}
	// Step 1: Classify characters into S-type (true) and L-type (false)
	t := make([]bool, n)
	t[n-1] = true // sentinel is S-type
	for i := n - 2; i >= 0; i-- {
		if s[i] < s[i+1] {
			t[i] = true
		} else if s[i] > s[i+1] {
			t[i] = false
		} else {
			t[i] = t[i+1]
		}
	}
	// Step 2: Identify LMS positions (positions i where t[i] is S-type and t[i-1] is L-type)
	var lmsPositions []int
	for i := 1; i < n; i++ {
		if t[i] && !t[i-1] {
			lmsPositions = append(lmsPositions, i)
		}
	}
	// Step 3: Induced sort the LMS positions
	SA = induceSort(s, SA, t, K, lmsPositions)
	// Step 4: Extract sorted LMS positions from SA
	var sortedLMS []int
	for _, pos := range SA {
		if pos > 0 && t[pos] && !t[pos-1] {
			sortedLMS = append(sortedLMS, pos)
		}
	}
	// Step 5: Name the LMS substrings
	lmsNames = lmsNames[:n]
	for i := range lmsNames {
		lmsNames[i] = -1
	}
	name := 0
	prev := -1
	for _, pos := range sortedLMS {
		if prev == -1 {
			lmsNames[pos] = name
		} else {
			if !lmsSubstringEqual(s, t, prev, pos) {
				name++
			}
			lmsNames[pos] = name
		}
		prev = pos
	}
	numNames := name + 1

	// Step 6: Build the reduced problem
	reduced := make([]int, 0, len(lmsPositions))
	for _, pos := range lmsPositions {
		reduced = append(reduced, lmsNames[pos])
	}

	var reducedSA []int
	if numNames < len(reduced) {
		reducedSA = SAIS(reduced, numNames, len(reduced), SA, lmsNames)
	} else {
		// If all names are unique then the order is already determined.
		reducedSA = make([]int, len(reduced))
		for i, name := range reduced {
			reducedSA[name] = i
		}
	}
	// Step 7: Map the reduced SA back to the original LMS positions and induce sort again.
	orderedLMS := make([]int, len(reducedSA))
	for i, idx := range reducedSA {
		orderedLMS[i] = lmsPositions[idx]
	}
	for i := range SA {
		SA[i] = -1
	}
	SA = induceSort(s, SA, t, K, orderedLMS)
	return SA
}

// induceSort performs the induced sorting given the LMS positions.
func induceSort(s []int, SA []int, t []bool, K int, lms []int) []int {
	bs := computeBucketSizes(s, K)
	bucketTails := computeBucketTails(bs)
	// Place LMS positions into the ends of their buckets.
	for i := len(lms) - 1; i >= 0; i-- {
		pos := lms[i]
		c := s[pos]
		SA[bucketTails[c]] = pos
		bucketTails[c]--
	}
	bucketHeads := computeBucketHeads(bs)
	// Induce L-type suffixes.
	for i := range SA {
		pos := SA[i]
		if pos > 0 && !t[pos-1] {
			c := s[pos-1]
			SA[bucketHeads[c]] = pos - 1
			bucketHeads[c]++
		}
	}
	bucketTails = computeBucketTails(bs)
	// Induce S-type suffixes.
	for i := len(SA) - 1; i >= 0; i-- {
		pos := SA[i]
		if pos > 0 && t[pos-1] {
			c := s[pos-1]
			SA[bucketTails[c]] = pos - 1
			bucketTails[c]--
		}
	}
	return SA
}

// computeBucketSizes returns a slice with the counts for each symbol.
func computeBucketSizes(s []int, K int) []int {
	bs := make([]int, K)
	for i := 0; i < len(s); i++ {
		bs[s[i]]++
	}
	return bs
}

// computeBucketHeads returns the starting index for each bucket.
func computeBucketHeads(bs []int) []int {
	heads := make([]int, len(bs))
	sum := 0
	for i, v := range bs {
		heads[i] = sum
		sum += v
	}
	return heads
}

// computeBucketTails returns the ending index for each bucket.
func computeBucketTails(bs []int) []int {
	tails := make([]int, len(bs))
	sum := 0
	for i, v := range bs {
		sum += v
		tails[i] = sum - 1
	}
	return tails
}

// lmsSubstringEqual compares the LMS substring starting at i with that at j.
func lmsSubstringEqual(s []int, t []bool, i, j int) bool {
	n := len(s)
	for {
		if s[i] != s[j] {
			return false
		}
		// Check if both positions are LMS positions.
		iIsLMS := (i > 0 && t[i] && !t[i-1])
		jIsLMS := (j > 0 && t[j] && !t[j-1])
		if iIsLMS && jIsLMS {
			return true
		}
		if iIsLMS != jIsLMS {
			return false
		}
		i++
		j++
		if i >= n || j >= n {
			break
		}
	}
	return false
}

// searchSequence uses binary search on the suffix array to locate all occurrences of query.
func searchSequence(genome string, sa []int, query string) []int {
	lb := lowerBound(genome, sa, query)
	if lb == -1 {
		return []int{}
	}
	ub := upperBound(genome, sa, query)
	return sa[lb:ub]
}

// lowerBound finds the first position in SA whose suffix is not less than query.
func lowerBound(genome string, sa []int, query string) int {
	lo := 0
	hi := len(sa)
	for lo < hi {
		mid := (lo + hi) / 2
		if compareSuffix(genome, sa[mid], query) < 0 {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	if lo < len(sa) && strings.HasPrefix(genome[sa[lo]:], query) {
		return lo
	}
	return -1
}

// upperBound finds the first position in SA whose suffix is greater than query.
func upperBound(genome string, sa []int, query string) int {
	lo := 0
	hi := len(sa)
	for lo < hi {
		mid := (lo + hi) / 2
		if compareSuffix(genome, sa[mid], query) <= 0 {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	return lo
}

// compareSuffix compares the suffix of genome starting at pos with query.
func compareSuffix(genome string, pos int, query string) int {
	i := 0
	for i < len(query) && pos+i < len(genome) {
		if genome[pos+i] != query[i] {
			return int(genome[pos+i]) - int(query[i])
		}
		i++
	}
	if i < len(query) {
		return -1
	}
	return 0
}
