var barrierList = [
{ name: "Short line",
    locations: [
12, 15,
12, 16,
12, 17,
12, 18,
12, 19,
12, 20,
12, 21,
12, 22,
12, 23]
},
{ name: "Long line",
    locations: [
13, 11,
13, 12,
13, 13,
13, 14,
13, 15,
13, 16,
13, 17,
13, 18,
13, 19,
13, 20,
13, 21,
13, 22,
13, 23,
13, 24,
13, 25,
13, 26,
13, 27,
13, 28
]
},
{ name: "Diagonal",
    locations: [
30, 14,
29, 15,
30, 15,
28, 16,
29, 16,
27, 17,
28, 17,
26, 18,
27, 18,
25, 19,
26, 19,
24, 20,
25, 20,
23, 21,
24, 21,
22, 22,
23, 22,
21, 23,
22, 23,
20, 24,
21, 24,
19, 25,
20, 25,
18, 26,
19, 26,
17, 27,
18, 27,
16, 28,
17, 28,
15, 29,
16, 29,
14, 30,
15, 30,
13, 31,
14, 31
]
},
{ name: "Shallow diagonal",
    locations: [
47, 18,
48, 18,
49, 18,
50, 18,
44, 19,
45, 19,
46, 19,
47, 19,
41, 20,
42, 20,
43, 20,
44, 20,
38, 21,
39, 21,
40, 21,
41, 21,
35, 22,
36, 22,
37, 22,
38, 22,
32, 23,
33, 23,
34, 23,
35, 23,
29, 24,
30, 24,
31, 24,
32, 24,
26, 25,
27, 25,
28, 25,
29, 25,
23, 26,
24, 26,
25, 26,
26, 26,
20, 27,
21, 27,
22, 27,
23, 27,
17, 28,
18, 28,
19, 28,
20, 28
]
},
{ name: "Small circle",
    locations: [
14, 11,
15, 11,
16, 11,
17, 11,
18, 11,
13, 12,
14, 12,
18, 12,
19, 12,
12, 13,
13, 13,
19, 13,
20, 13,
12, 14,
20, 14,
12, 15,
20, 15,
12, 16,
20, 16,
12, 17,
13, 17,
19, 17,
20, 17,
13, 18,
14, 18,
18, 18,
19, 18,
14, 19,
15, 19,
16, 19,
17, 19,
18, 19
]
},
{ name: "Large circle",
    locations: [
19, 11,
20, 11,
21, 11,
22, 11,
23, 11,
24, 11,
17, 12,
18, 12,
19, 12,
24, 12,
25, 12,
26, 12,
16, 13,
17, 13,
26, 13,
27, 13,
15, 14,
16, 14,
27, 14,
28, 14,
14, 15,
15, 15,
28, 15,
29, 15,
14, 16,
29, 16,
13, 17,
14, 17,
29, 17,
30, 17,
13, 18,
30, 18,
13, 19,
30, 19,
13, 20,
30, 20,
13, 21,
30, 21,
13, 22,
14, 22,
29, 22,
30, 22,
14, 23,
29, 23,
14, 24,
15, 24,
28, 24,
29, 24,
15, 25,
16, 25,
27, 25,
28, 25,
16, 26,
17, 26,
26, 26,
27, 26,
17, 27,
18, 27,
19, 27,
24, 27,
25, 27,
26, 27,
19, 28,
20, 28,
21, 28,
22, 28,
23, 28,
24, 28
]
},
{ name: "Line with spoiler",
    locations: [
16, 20,
16, 21,
16, 22,
16, 23,
16, 24,
17, 24,
18, 24,
19, 24,
20, 24,
21, 24,
22, 24,
23, 24,
24, 24,
25, 24,
26, 24,
27, 24,
28, 24,
29, 24,
30, 24,
31, 24,
32, 24,
33, 24,
34, 24,
35, 24,
36, 24,
37, 24,
38, 24,
39, 24,
40, 24,
41, 24,
42, 24,
43, 24,
44, 24,
45, 24,
46, 24,
47, 24,
48, 24,
49, 24,
50, 24,
16, 25,
16, 26,
16, 27,
16, 28
]
},
{ name: "Circle with spoiler",
    locations: [
29, 36,
30, 36,
31, 36,
32, 36,
33, 36,
28, 37,
29, 37,
33, 37,
34, 37,
27, 38,
28, 38,
34, 38,
35, 38,
27, 39,
35, 39,
27, 40,
35, 40,
36, 40,
37, 40,
38, 40,
39, 40,
40, 40,
41, 40,
42, 40,
43, 40,
44, 40,
45, 40,
46, 40,
47, 40,
48, 40,
49, 40,
50, 40,
51, 40,
52, 40,
53, 40,
54, 40,
55, 40,
56, 40,
57, 40,
58, 40,
59, 40,
60, 40,
61, 40,
62, 40,
63, 40,
64, 40,
65, 40,
66, 40,
67, 40,
68, 40,
69, 40,
27, 41,
35, 41,
27, 42,
28, 42,
34, 42,
35, 42,
28, 43,
29, 43,
33, 43,
34, 43,
29, 44,
30, 44,
31, 44,
32, 44,
33, 44
]
},
{ name: "Right angle",
    locations: [
27, 36,
28, 36,
29, 36,
30, 36,
31, 36,
32, 36,
33, 36,
34, 36,
35, 36,
36, 36,
37, 36,
38, 36,
39, 36,
40, 36,
41, 36,
42, 36,
43, 36,
44, 36,
45, 36,
46, 36,
47, 36,
48, 36,
49, 36,
50, 36,
51, 36,
52, 36,
53, 36,
54, 36,
55, 36,
56, 36,
57, 36,
58, 36,
59, 36,
60, 36,
61, 36,
62, 36,
63, 36,
64, 36,
65, 36,
66, 36,
67, 36,
68, 36,
69, 36,
70, 36,
71, 36,
72, 36,
73, 36,
74, 36,
75, 36,
76, 36,
77, 36,
78, 36,
79, 36,
27, 37,
27, 38,
27, 39,
27, 40,
27, 41,
27, 42,
27, 43,
27, 44
]
},
{ name: "Wedge",
    locations: [
27, 36,
28, 36,
29, 36,
30, 36,
31, 36,
32, 36,
33, 36,
34, 36,
35, 36,
36, 36,
37, 36,
38, 36,
39, 36,
40, 36,
41, 36,
42, 36,
43, 36,
44, 36,
45, 36,
46, 36,
47, 36,
48, 36,
49, 36,
50, 36,
51, 36,
52, 36,
53, 36,
54, 36,
55, 36,
56, 36,
57, 36,
58, 36,
59, 36,
60, 36,
61, 36,
62, 36,
63, 36,
64, 36,
65, 36,
66, 36,
67, 36,
68, 36,
69, 36,
70, 36,
71, 36,
72, 36,
73, 36,
74, 36,
75, 36,
76, 36,
77, 36,
78, 36,
79, 36,
27, 37,
67, 37,
68, 37,
69, 37,
70, 37,
71, 37,
72, 37,
73, 37,
27, 38,
61, 38,
62, 38,
63, 38,
64, 38,
65, 38,
66, 38,
67, 38,
27, 39,
55, 39,
56, 39,
57, 39,
58, 39,
59, 39,
60, 39,
61, 39,
27, 40,
49, 40,
50, 40,
51, 40,
52, 40,
53, 40,
54, 40,
55, 40,
27, 41,
43, 41,
44, 41,
45, 41,
46, 41,
47, 41,
48, 41,
49, 41,
27, 42,
37, 42,
38, 42,
39, 42,
40, 42,
41, 42,
42, 42,
43, 42,
27, 43,
31, 43,
32, 43,
33, 43,
34, 43,
35, 43,
36, 43,
37, 43,
27, 44,
28, 44,
29, 44,
30, 44,
31, 44
]
},
{ name: "Airfoil",
    locations: [
17, 16,
18, 16,
19, 16,
20, 16,
21, 16,
22, 16,
23, 16,
24, 16,
25, 16,
26, 16,
27, 16,
28, 16,
29, 16,
30, 16,
31, 16,
32, 16,
33, 16,
34, 16,
35, 16,
36, 16,
37, 16,
38, 16,
39, 16,
40, 16,
41, 16,
42, 16,
43, 16,
44, 16,
45, 16,
46, 16,
47, 16,
48, 16,
49, 16,
50, 16,
51, 16,
52, 16,
53, 16,
54, 16,
55, 16,
56, 16,
57, 16,
58, 16,
59, 16,
60, 16,
61, 16,
62, 16,
63, 16,
64, 16,
65, 16,
66, 16,
67, 16,
68, 16,
14, 17,
15, 17,
16, 17,
17, 17,
56, 17,
57, 17,
58, 17,
59, 17,
60, 17,
61, 17,
62, 17,
13, 18,
14, 18,
50, 18,
51, 18,
52, 18,
53, 18,
54, 18,
55, 18,
56, 18,
13, 19,
44, 19,
45, 19,
46, 19,
47, 19,
48, 19,
49, 19,
50, 19,
13, 20,
38, 20,
39, 20,
40, 20,
41, 20,
42, 20,
43, 20,
44, 20,
13, 21,
14, 21,
32, 21,
33, 21,
34, 21,
35, 21,
36, 21,
37, 21,
38, 21,
14, 22,
15, 22,
26, 22,
27, 22,
28, 22,
29, 22,
30, 22,
31, 22,
32, 22,
15, 23,
16, 23,
17, 23,
18, 23,
21, 23,
22, 23,
23, 23,
24, 23,
25, 23,
26, 23,
18, 24,
19, 24,
20, 24,
21, 24
]
}
];