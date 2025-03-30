package main

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
