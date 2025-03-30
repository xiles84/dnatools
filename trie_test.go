package main

import (
	"reflect"
	"testing"
)

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
