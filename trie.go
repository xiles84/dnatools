package main

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
