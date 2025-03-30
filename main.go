package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// SuffixEntry holds the suffix array entry, the originating line, and LCP value.
type SuffixEntry struct {
	Pos  int
	Line int
	LCP  int
}

// TrieNode represents a node in the trie.
type TrieNode struct {
	children map[rune]*TrieNode
	isEnd    bool
	pattern  string
}

// NewTrie creates a new empty trie node.
func NewTrie() *TrieNode {
	return &TrieNode{children: make(map[rune]*TrieNode)}
}

// Insert adds a pattern into the trie.
func (node *TrieNode) Insert(pattern string) {
	current := node
	for _, ch := range pattern {
		if current.children == nil {
			current.children = make(map[rune]*TrieNode)
		}
		if _, exists := current.children[ch]; !exists {
			current.children[ch] = &TrieNode{children: make(map[rune]*TrieNode)}
		}
		current = current.children[ch]
	}
	current.isEnd = true
	current.pattern = pattern
}

// searchTrie scans the text and returns a map where each key is a pattern found
// and the value is a slice of starting positions where that pattern occurs.
func searchTrie(text string, root *TrieNode) map[string][]int {
	results := make(map[string][]int)
	// For each position in the text, try to match a pattern.
	for i := 0; i < len(text); i++ {
		current := root
		for j := i; j < len(text); j++ {
			ch := rune(text[j])
			next, exists := current.children[ch]
			if !exists {
				break
			}
			current = next
			if current.isEnd {
				results[current.pattern] = append(results[current.pattern], i)
			}
		}
	}
	return results
}

func main() {
	indexMode := flag.Bool("m", false, "Index mode: build suffix array index")
	searchQueryStr := flag.String("s", "", "Search mode: search for sequence using suffix array")
	trieFile := flag.String("t", "", "Trie search mode: file containing multiple query patterns (one per line)")
	fileName := flag.String("f", "genoma.txt", "Genome file name")
	flag.Parse()

	// Read genome file.
	data, err := os.ReadFile(*fileName)
	if err != nil {
		fmt.Println("Error reading genome file:", err)
		os.Exit(1)
	}

	// Build genome from file lines. Each nonempty line is a DNA sequence.
	lines := strings.Split(string(data), "\n")
	var sequences []string
	for _, line := range lines {
		trimmed := strings.TrimSpace(line)
		if trimmed != "" {
			sequences = append(sequences, trimmed)
		}
	}
	var genomeBuilder strings.Builder
	var lineMap []int
	for i, seq := range sequences {
		for _, ch := range seq {
			genomeBuilder.WriteRune(ch)
			lineMap = append(lineMap, i) // record originating line index
		}
		// Append a separator if not the last sequence.
		if i < len(sequences)-1 {
			genomeBuilder.WriteByte('$')
			lineMap = append(lineMap, -1) // -1 indicates separator
		}
	}
	genome := genomeBuilder.String()

	// Index mode using suffix array (with LCP) remains as before.
	if *indexMode {
		fmt.Println("Building suffix array index using SAIS algorithm...")
		encoded, alphabetSize := encodeString(genome)
		sa := SAISEntryPoint(encoded, alphabetSize)
		lcp := computeLCP(genome, sa)
		entries := make([]SuffixEntry, len(sa))
		for i, pos := range sa {
			lineNum := -1
			if pos < len(lineMap) {
				lineNum = lineMap[pos]
			}
			entries[i] = SuffixEntry{Pos: pos, Line: lineNum, LCP: lcp[i]}
		}
		err = saveIndex("sa.idx", entries)
		if err != nil {
			fmt.Println("Error saving index:", err)
			os.Exit(1)
		}
		fmt.Println("Index built and saved to sa.idx")
		// Suffix array search mode.
	} else if *searchQueryStr != "" {
		entries, err := loadIndex("sa.idx")
		if err != nil {
			fmt.Println("Error loading index:", err)
			os.Exit(1)
		}
		fmt.Printf("Searching for sequence: %s\n", *searchQueryStr)
		results := searchSequence(genome, entries, *searchQueryStr)
		if len(results) == 0 {
			fmt.Println("Sequence not found.")
		} else {
			fmt.Println("Sequence found at positions (global position, DNA line):")
			for _, entry := range results {
				fmt.Printf("(%d, %d) ", entry.Pos, entry.Line)
			}
			fmt.Println()
		}
		// Trie search mode: used with the -t flag.
	} else if *trieFile != "" {
		// Read the file containing multiple query patterns.
		patternData, err := os.ReadFile(*trieFile)
		if err != nil {
			fmt.Println("Error reading trie file:", err)
			os.Exit(1)
		}
		patternLines := strings.Split(string(patternData), "\n")
		var patterns []string
		for _, l := range patternLines {
			trimmed := strings.TrimSpace(l)
			if trimmed != "" {
				patterns = append(patterns, trimmed)
			}
		}
		// Build the trie from the query patterns.
		trie := NewTrie()
		for _, pat := range patterns {
			trie.Insert(pat)
		}
		// Search the genome using the trie.
		results := searchTrie(genome, trie)
		// Print results.
		for pat, positions := range results {
			// Optionally, you could annotate the positions with DNA line numbers from lineMap.
			fmt.Printf("Pattern %q found at positions: %v\n", pat, positions)
		}
	} else {
		fmt.Println("Please provide -m to build index, -s <sequence> for suffix array search, or -t <file> for trie search (used with -f).")
	}
}

// saveIndex writes the suffix entries to a file (one entry per line: pos, line, LCP).
func saveIndex(filename string, entries []SuffixEntry) error {
	var lines []string
	for _, entry := range entries {
		line := fmt.Sprintf("%d %d %d", entry.Pos, entry.Line, entry.LCP)
		lines = append(lines, line)
	}
	content := strings.Join(lines, "\n")
	return os.WriteFile(filename, []byte(content), 0644)
}

// loadIndex reads the suffix entries from a file.
func loadIndex(filename string) ([]SuffixEntry, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	var entries []SuffixEntry
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		parts := strings.Fields(line)
		if len(parts) != 3 {
			continue
		}
		pos, err1 := strconv.Atoi(parts[0])
		lineNum, err2 := strconv.Atoi(parts[1])
		lcpVal, err3 := strconv.Atoi(parts[2])
		if err1 != nil || err2 != nil || err3 != nil {
			continue
		}
		entries = append(entries, SuffixEntry{Pos: pos, Line: lineNum, LCP: lcpVal})
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return entries, nil
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
	t := make([]bool, n)
	t[n-1] = true
	for i := n - 2; i >= 0; i-- {
		if s[i] < s[i+1] {
			t[i] = true
		} else if s[i] > s[i+1] {
			t[i] = false
		} else {
			t[i] = t[i+1]
		}
	}
	var lmsPositions []int
	for i := 1; i < n; i++ {
		if t[i] && !t[i-1] {
			lmsPositions = append(lmsPositions, i)
		}
	}
	SA = induceSort(s, SA, t, K, lmsPositions)
	var sortedLMS []int
	for _, pos := range SA {
		if pos > 0 && t[pos] && !t[pos-1] {
			sortedLMS = append(sortedLMS, pos)
		}
	}
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
	reduced := make([]int, 0, len(lmsPositions))
	for _, pos := range lmsPositions {
		reduced = append(reduced, lmsNames[pos])
	}
	var reducedSA []int
	if numNames < len(reduced) {
		reducedSA = SAIS(reduced, numNames, len(reduced), SA, lmsNames)
	} else {
		reducedSA = make([]int, len(reduced))
		for i, name := range reduced {
			reducedSA[name] = i
		}
	}
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

func induceSort(s []int, SA []int, t []bool, K int, lms []int) []int {
	bs := computeBucketSizes(s, K)
	bucketTails := computeBucketTails(bs)
	for i := len(lms) - 1; i >= 0; i-- {
		pos := lms[i]
		c := s[pos]
		SA[bucketTails[c]] = pos
		bucketTails[c]--
	}
	bucketHeads := computeBucketHeads(bs)
	for i := range SA {
		pos := SA[i]
		if pos > 0 && !t[pos-1] {
			c := s[pos-1]
			SA[bucketHeads[c]] = pos - 1
			bucketHeads[c]++
		}
	}
	bucketTails = computeBucketTails(bs)
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

func computeBucketSizes(s []int, K int) []int {
	bs := make([]int, K)
	for i := 0; i < len(s); i++ {
		bs[s[i]]++
	}
	return bs
}

func computeBucketHeads(bs []int) []int {
	heads := make([]int, len(bs))
	sum := 0
	for i, v := range bs {
		heads[i] = sum
		sum += v
	}
	return heads
}

func computeBucketTails(bs []int) []int {
	tails := make([]int, len(bs))
	sum := 0
	for i, v := range bs {
		sum += v
		tails[i] = sum - 1
	}
	return tails
}

func lmsSubstringEqual(s []int, t []bool, i, j int) bool {
	n := len(s)
	for {
		if s[i] != s[j] {
			return false
		}
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

// computeLCP computes the Longest Common Prefix array using Kasai's algorithm.
func computeLCP(s string, sa []int) []int {
	n := len(sa)
	lcp := make([]int, n)
	rank := make([]int, n)
	for i, pos := range sa {
		if pos < n {
			rank[pos] = i
		}
	}
	h := 0
	for i := 0; i < n; i++ {
		if rank[i] > 0 {
			j := sa[rank[i]-1]
			for i+h < len(s) && j+h < len(s) && s[i+h] == s[j+h] {
				h++
			}
			lcp[rank[i]] = h
			if h > 0 {
				h--
			}
		} else {
			lcp[rank[i]] = 0
		}
	}
	return lcp
}

// searchSequence uses binary search on the suffix entries to locate all occurrences of query.
func searchSequence(genome string, entries []SuffixEntry, query string) []SuffixEntry {
	lb := lowerBound(genome, entries, query)
	if lb == -1 {
		return []SuffixEntry{}
	}
	ub := upperBound(genome, entries, query)
	return entries[lb:ub]
}

func lowerBound(genome string, entries []SuffixEntry, query string) int {
	lo := 0
	hi := len(entries)
	for lo < hi {
		mid := (lo + hi) / 2
		if compareSuffix(genome, entries[mid].Pos, query) < 0 {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	if lo < len(entries) && strings.HasPrefix(genome[entries[lo].Pos:], query) {
		return lo
	}
	return -1
}

func upperBound(genome string, entries []SuffixEntry, query string) int {
	lo := 0
	hi := len(entries)
	for lo < hi {
		mid := (lo + hi) / 2
		if compareSuffix(genome, entries[mid].Pos, query) <= 0 {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	return lo
}

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
