package main

import (
	"os"
	"reflect"
	"strings"
	"testing"
)

func TestSuffixArray(t *testing.T) {
	// Test cases for SAISEntryPoint on single strings.
	testCases := []struct {
		input    string
		expected []int
	}{
		{
			input:    "banana",
			expected: []int{6, 5, 3, 1, 0, 4, 2},
		},
		{
			input:    "mississippi",
			expected: []int{11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2},
		},
		{
			input:    "a",
			expected: []int{1, 0},
		},
		{
			input:    "",
			expected: []int{0},
		},
	}

	for _, tc := range testCases {
		encoded, alphabetSize := encodeString(tc.input)
		sa := SAISEntryPoint(encoded, alphabetSize)
		if !reflect.DeepEqual(sa, tc.expected) {
			t.Errorf("Suffix array for %q: got %v, expected %v", tc.input, sa, tc.expected)
		}
	}
}

func TestComputeLCP(t *testing.T) {
	// For "banana", the expected suffix array is [6,5,3,1,0,4,2]
	// and the expected LCP array is [0, 0, 1, 3, 0, 0, 2]
	input := "banana"
	encoded, alphabetSize := encodeString(input)
	sa := SAISEntryPoint(encoded, alphabetSize)
	// Note: We pass input + "$" to simulate the trailing sentinel.
	lcp := computeLCP(input+"$", sa)
	expectedLCP := []int{0, 0, 1, 3, 0, 0, 2}
	if !reflect.DeepEqual(lcp, expectedLCP) {
		t.Errorf("LCP array for %q: got %v, expected %v", input, lcp, expectedLCP)
	}
}

func TestEnhancedIndex(t *testing.T) {
	// Simulate a file with two DNA sequences.
	sequences := []string{"ACGT", "TGCA"}
	var genomeBuilder strings.Builder
	var lineMap []int
	for i, seq := range sequences {
		for _, ch := range seq {
			genomeBuilder.WriteRune(ch)
			lineMap = append(lineMap, i) // record originating sequence index
		}
		if i < len(sequences)-1 {
			genomeBuilder.WriteByte('$')  // separator character
			lineMap = append(lineMap, -1) // -1 indicates separator
		}
	}
	genome := genomeBuilder.String()

	encoded, alphabetSize := encodeString(genome)
	sa := SAISEntryPoint(encoded, alphabetSize)
	lcp := computeLCP(genome, sa)

	// Build suffix entries containing global position, originating DNA line, and LCP.
	entries := make([]SuffixEntry, len(sa))
	for i, pos := range sa {
		lineNum := -1
		if pos < len(lineMap) {
			lineNum = lineMap[pos]
		}
		entries[i] = SuffixEntry{Pos: pos, Line: lineNum, LCP: lcp[i]}
	}

	// Verify that each entry’s line mapping is correct.
	for _, entry := range entries {
		if entry.Pos < len(lineMap) && entry.Line != lineMap[entry.Pos] {
			t.Errorf("For global position %d, expected line %d, got %d", entry.Pos, lineMap[entry.Pos], entry.Line)
		}
	}

	// Additionally, test the save/load index functionality.
	tmpFile := "test_sa.idx"
	defer os.Remove(tmpFile)
	if err := saveIndex(tmpFile, entries); err != nil {
		t.Fatalf("Failed to save index: %v", err)
	}
	loadedEntries, err := loadIndex(tmpFile)
	if err != nil {
		t.Fatalf("Failed to load index: %v", err)
	}
	if !reflect.DeepEqual(loadedEntries, entries) {
		t.Errorf("Loaded entries do not match. Got %v, expected %v", loadedEntries, entries)
	}
}

func TestSearchSequence(t *testing.T) {
	// Build an index from two DNA sequences.
	sequences := []string{"ACGT", "TGCA"}
	var genomeBuilder strings.Builder
	var lineMap []int
	for i, seq := range sequences {
		for _, ch := range seq {
			genomeBuilder.WriteRune(ch)
			lineMap = append(lineMap, i)
		}
		if i < len(sequences)-1 {
			genomeBuilder.WriteByte('$')
			lineMap = append(lineMap, -1)
		}
	}
	genome := genomeBuilder.String()

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

	// Search for a query that should occur in the first sequence ("ACGT").
	query := "CG"
	results := searchSequence(genome, entries, query)
	if len(results) == 0 {
		t.Errorf("Expected to find query %q in genome", query)
	}
	for _, entry := range results {
		if entry.Line != 0 {
			t.Errorf("Query %q expected in line 0 but found in line %d", query, entry.Line)
		}
	}

	// Search for a query that should occur in the second sequence ("TGCA").
	query2 := "GC"
	results2 := searchSequence(genome, entries, query2)
	if len(results2) == 0 {
		t.Errorf("Expected to find query %q in genome", query2)
	}
	for _, entry := range results2 {
		if entry.Line != 1 {
			t.Errorf("Query %q expected in line 1 but found in line %d", query2, entry.Line)
		}
	}
}

func TestTrieSearch(t *testing.T) {
	// Sample genome string to search.
	text := "ACGTACGT"
	// Build a trie with multiple patterns.
	trie := NewTrie()
	patterns := []string{"ACG", "CGT", "TAC", "GTAC"}
	for _, pat := range patterns {
		trie.Insert(pat)
	}

	// Execute the trie search on the genome.
	results := searchTrie(text, trie)
	// Expected occurrences:
	// "ACG" should be found starting at indices 0 and 4.
	// "CGT" should be found starting at indices 1 and 5.
	// "TAC" should be found starting at index 3.
	// "GTAC" should be found starting at index 2.
	expected := map[string][]int{
		"ACG":  {0, 4},
		"CGT":  {1, 5},
		"TAC":  {3},
		"GTAC": {2},
	}
	if !reflect.DeepEqual(results, expected) {
		t.Errorf("Trie search results mismatch. Expected %v, got %v", expected, results)
	}
}

func BenchmarkSAIS(b *testing.B) {
	// Generate a synthetic genome sequence for benchmarking.
	genome := strings.Repeat("ACGT", 10000000) // 40,000,000 characters
	encoded, alphabetSize := encodeString(genome)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = SAISEntryPoint(encoded, alphabetSize)
	}
}
