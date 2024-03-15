/*
 *     Author: Philip R. Kensche, DKFZ Heidelberg (Omics IT and Data Management Core Facility)
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *     LICENSE: GPL
 */

#include <vector>
#include <string>
#include <stdexcept>
#include <boost/exception/all.hpp>

#include "Hg37ChrConverter.h"
#include "global.h"
#include "IndexRange.h"


namespace sophia {

    namespace hg37 {

        static const std::vector<ChrName> indexToChrName {
            "0",          "1",          "2",          "3",          "4",
            "5",          "6",          "7",          "8",          "9",
            "10",         "11",         "12",         "13",         "14",
            "15",         "16",         "17",         "18",         "19",
            "20",         "21",         "22",         "23",         "24",
            "25",         "26",         "27",         "28",         "29",
            "30",         "31",         "32",         "33",         "34",
            "35",         "36",         "37",         "38",         "39",
            "X",          "Y",          "42",         "43",         "44",
            "45",         "46",         "47",         "48",         "49",
            "50",         "51",         "52",         "53",         "54",
            "55",         "56",         "57",         "58",         "59",
            "60",         "61",         "62",         "63",         "64",
            "65",         "66",         "67",         "68",         "69",
            "70",         "71",         "72",         "73",         "74",
            "75",         "76",         "77",         "78",         "79",
            "80",         "81",         "82",         "83",         "84",
            "85",         "86",         "87",         "88",         "89",
            "90",         "91",         "92",         "93",         "94",
            "95",         "96",         "97",         "98",         "99",
            "100",        "101",        "102",        "103",        "104",
            "105",        "106",        "107",        "108",        "109",
            "110",        "111",        "112",        "113",        "114",
            "115",        "116",        "117",        "118",        "119",
            "120",        "121",        "122",        "123",        "124",
            "125",        "126",        "127",        "128",        "129",
            "130",        "131",        "132",        "133",        "134",
            "135",        "136",        "137",        "138",        "139",
            "140",        "141",        "142",        "143",        "144",
            "145",        "146",        "147",        "148",        "149",
            "150",        "151",        "152",        "153",        "154",
            "155",        "156",        "157",        "158",        "159",
            "160",        "161",        "162",        "163",        "164",
            "165",        "166",        "167",        "168",        "169",
            "170",        "171",        "172",        "173",        "174",
            "175",        "176",        "177",        "178",        "179",
            "180",        "181",        "182",        "183",        "184",
            "185",        "186",        "187",        "188",        "189",
            "190",        "GL000191.1", "GL000192.1", "GL000193.1", "GL000194.1",
            "GL000195.1", "GL000196.1", "GL000197.1", "GL000198.1", "GL000199.1",
            "GL000200.1", "GL000201.1", "GL000202.1", "GL000203.1", "GL000204.1",
            "GL000205.1", "GL000206.1", "GL000207.1", "GL000208.1", "GL000209.1",
            "GL000210.1", "GL000211.1", "GL000212.1", "GL000213.1", "GL000214.1",
            "GL000215.1", "GL000216.1", "GL000217.1", "GL000218.1", "GL000219.1",
            "GL000220.1", "GL000221.1", "GL000222.1", "GL000223.1", "GL000224.1",
            "GL000225.1", "GL000226.1", "GL000227.1", "GL000228.1", "GL000229.1",
            "GL000230.1", "GL000231.1", "GL000232.1", "GL000233.1", "GL000234.1",
            "GL000235.1", "GL000236.1", "GL000237.1", "GL000238.1", "GL000239.1",
            "GL000240.1", "GL000241.1", "GL000242.1", "GL000243.1", "GL000244.1",
            "GL000245.1", "GL000246.1", "GL000247.1", "GL000248.1", "GL000249.1",
            "250",        "251",        "252",        "253",        "254",
            "255",        "256",        "257",        "258",        "259",
            "260",        "261",        "262",        "263",        "264",
            "265",        "266",        "267",        "268",        "269",
            "270",        "271",        "272",        "273",        "274",
            "275",        "276",        "277",        "278",        "279",
            "280",        "281",        "282",        "283",        "284",
            "285",        "286",        "287",        "288",        "289",
            "290",        "291",        "292",        "293",        "294",
            "295",        "296",        "297",        "298",        "299",
            "300",        "301",        "302",        "303",        "304",
            "305",        "306",        "307",        "308",        "309",
            "310",        "311",        "312",        "313",        "314",
            "315",        "316",        "317",        "318",        "319",
            "320",        "321",        "322",        "323",        "324",
            "325",        "326",        "327",        "328",        "329",
            "330",        "331",        "332",        "333",        "334",
            "335",        "336",        "337",        "338",        "339",
            "340",        "341",        "342",        "343",        "344",
            "345",        "346",        "347",        "348",        "349",
            "350",        "351",        "352",        "353",        "354",
            "355",        "356",        "357",        "358",        "359",
            "360",        "361",        "362",        "363",        "364",
            "365",        "366",        "367",        "368",        "369",
            "370",        "371",        "372",        "373",        "374",
            "375",        "376",        "377",        "378",        "379",
            "380",        "381",        "382",        "383",        "384",
            "385",        "386",        "387",        "388",        "389",
            "390",        "391",        "392",        "393",        "394",
            "395",        "396",        "397",        "398",        "399",
            "400",        "401",        "402",        "403",        "404",
            "405",        "406",        "407",        "408",        "409",
            "410",        "411",        "412",        "413",        "414",
            "415",        "416",        "417",        "418",        "419",
            "420",        "421",        "422",        "423",        "424",
            "425",        "426",        "427",        "428",        "429",
            "430",        "431",        "432",        "433",        "434",
            "435",        "436",        "437",        "438",        "439",
            "440",        "441",        "442",        "443",        "444",
            "445",        "446",        "447",        "448",        "449",
            "450",        "451",        "452",        "453",        "454",
            "455",        "456",        "457",        "458",        "459",
            "460",        "461",        "462",        "463",        "464",
            "465",        "466",        "467",        "468",        "469",
            "470",        "471",        "472",        "473",        "474",
            "475",        "476",        "477",        "478",        "479",
            "480",        "481",        "482",        "483",        "484",
            "485",        "486",        "487",        "488",        "489",
            "490",        "491",        "492",        "493",        "494",
            "495",        "496",        "497",        "498",        "499",
            "500",        "501",        "502",        "503",        "504",
            "505",        "506",        "507",        "508",        "509",
            "510",        "511",        "512",        "513",        "514",
            "515",        "516",        "517",        "518",        "519",
            "520",        "521",        "522",        "523",        "524",
            "525",        "526",        "527",        "528",        "529",
            "530",        "531",        "532",        "533",        "534",
            "535",        "536",        "537",        "538",        "539",
            "540",        "541",        "542",        "543",        "544",
            "545",        "546",        "547",        "548",        "549",
            "550",        "551",        "552",        "553",        "554",
            "555",        "556",        "557",        "558",        "559",
            "560",        "561",        "562",        "563",        "564",
            "565",        "566",        "567",        "568",        "569",
            "570",        "571",        "572",        "573",        "574",
            "575",        "576",        "577",        "578",        "579",
            "580",        "581",        "582",        "583",        "584",
            "585",        "586",        "587",        "588",        "589",
            "590",        "591",        "592",        "593",        "594",
            "595",        "596",        "597",        "598",        "599",
            "600",        "601",        "602",        "603",        "604",
            "605",        "606",        "607",        "608",        "609",
            "610",        "611",        "612",        "613",        "614",
            "615",        "616",        "617",        "618",        "619",
            "620",        "621",        "622",        "623",        "624",
            "625",        "626",        "627",        "628",        "629",
            "630",        "631",        "632",        "633",        "634",
            "635",        "636",        "637",        "638",        "639",
            "640",        "641",        "642",        "643",        "644",
            "645",        "646",        "647",        "648",        "649",
            "650",        "651",        "652",        "653",        "654",
            "655",        "656",        "657",        "658",        "659",
            "660",        "661",        "662",        "663",        "664",
            "665",        "666",        "667",        "668",        "669",
            "670",        "671",        "672",        "673",        "674",
            "675",        "676",        "677",        "678",        "679",
            "680",        "681",        "682",        "683",        "684",
            "685",        "686",        "687",        "688",        "689",
            "690",        "691",        "692",        "693",        "694",
            "695",        "696",        "697",        "698",        "699",
            "700",        "701",        "702",        "703",        "704",
            "705",        "706",        "707",        "708",        "709",
            "710",        "711",        "712",        "713",        "714",
            "715",        "716",        "717",        "718",        "719",
            "720",        "721",        "722",        "723",        "724",
            "725",        "726",        "727",        "728",        "729",
            "730",        "731",        "732",        "733",        "734",
            "735",        "736",        "737",        "738",        "739",
            "740",        "741",        "742",        "743",        "744",
            "745",        "746",        "747",        "748",        "749",
            "750",        "751",        "752",        "753",        "754",
            "755",        "756",        "757",        "758",        "759",
            "760",        "761",        "762",        "763",        "764",
            "765",        "766",        "767",        "768",        "769",
            "770",        "771",        "772",        "773",        "774",
            "775",        "776",        "777",        "778",        "779",
            "780",        "781",        "782",        "783",        "784",
            "785",        "786",        "787",        "788",        "789",
            "790",        "791",        "792",        "793",        "794",
            "795",        "796",        "797",        "798",        "799",
            "800",        "801",        "802",        "803",        "804",
            "805",        "806",        "807",        "808",        "809",
            "810",        "811",        "812",        "813",        "814",
            "815",        "816",        "817",        "818",        "819",
            "820",        "821",        "822",        "823",        "824",
            "825",        "826",        "827",        "828",        "829",
            "830",        "831",        "832",        "833",        "834",
            "835",        "836",        "837",        "838",        "839",
            "840",        "841",        "842",        "843",        "844",
            "845",        "846",        "847",        "848",        "849",
            "850",        "851",        "852",        "853",        "854",
            "855",        "856",        "857",        "858",        "859",
            "860",        "861",        "862",        "863",        "864",
            "865",        "866",        "867",        "868",        "869",
            "870",        "871",        "872",        "873",        "874",
            "875",        "876",        "877",        "878",        "879",
            "880",        "881",        "882",        "883",        "884",
            "885",        "886",        "887",        "888",        "889",
            "890",        "891",        "892",        "893",        "894",
            "895",        "896",        "897",        "898",        "899",
            "900",        "901",        "902",        "903",        "904",
            "905",        "906",        "907",        "908",        "909",
            "910",        "911",        "912",        "913",        "914",
            "915",        "916",        "917",        "918",        "919",
            "920",        "921",        "922",        "923",        "924",
            "925",        "926",        "927",        "928",        "929",
            "930",        "931",        "932",        "933",        "934",
            "935",        "936",        "937",        "938",        "939",
            "940",        "941",        "942",        "943",        "944",
            "945",        "946",        "947",        "948",        "949",
            "950",        "951",        "952",        "953",        "954",
            "955",        "956",        "957",        "958",        "959",
            "960",        "961",        "962",        "963",        "964",
            "965",        "966",        "967",        "968",        "969",
            "970",        "971",        "972",        "973",        "974",
            "975",        "976",        "977",        "978",        "979",
            "980",        "981",        "982",        "983",        "984",
            "985",        "986",        "987",        "988",        "989",
            "990",        "991",        "992",        "993",        "994",
            "995",        "996",        "997",        "998",        "hs37d5",
            "NC_007605",  "MT",         "phiX174",    "INVALID"};

        static const ChrIndex ZERO = 0;
        static const ChrIndex xIndex = 40;
        static const ChrIndex yIndex = 41;
        static const ChrIndex decoyIndex = 999;
        static const ChrIndex virusIndex = 1000;
        static const ChrIndex mtIndex = 1001;
        static const ChrIndex phixIndex = 1002;
        static const ChrIndex INVALID = 1003;

        static const IndexRange automoseRange = {1, 23};
        static const IndexRange gonosomeRange = {40, 42};
        static const IndexRange unassignedRange = {191, 250};
        static const IndexRange decoyRange = {decoyIndex, decoyIndex + 1};
        static const IndexRange virusRange = {virusIndex, virusIndex + 1};
        static const IndexRange extrachromosomalRange = {mtIndex, mtIndex + 1};
        static const IndexRange technicalRange = {phixIndex, phixIndex + 1};

        /* 85 compressed mref chromosomes */
        static const std::vector<ChrName> compressedMrefIndexToChrName {
            "1",          "2",          "3",          "4",          "5",
            "6",          "7",          "8",          "9",          "10",
            "11",         "12",         "13",         "14",         "15",
            "16",         "17",         "18",         "19",         "20",
            "21",         "22",         "X",          "Y",          "GL000191.1",
            "GL000192.1", "GL000193.1", "GL000194.1", "GL000195.1", "GL000196.1",
            "GL000197.1", "GL000198.1", "GL000199.1", "GL000200.1", "GL000201.1",
            "GL000202.1", "GL000203.1", "GL000204.1", "GL000205.1", "GL000206.1",
            "GL000207.1", "GL000208.1", "GL000209.1", "GL000210.1", "GL000211.1",
            "GL000212.1", "GL000213.1", "GL000214.1", "GL000215.1", "GL000216.1",
            "GL000217.1", "GL000218.1", "GL000219.1", "GL000220.1", "GL000221.1",
            "GL000222.1", "GL000223.1", "GL000224.1", "GL000225.1", "GL000226.1",
            "GL000227.1", "GL000228.1", "GL000229.1", "GL000230.1", "GL000231.1",
            "GL000232.1", "GL000233.1", "GL000234.1", "GL000235.1", "GL000236.1",
            "GL000237.1", "GL000238.1", "GL000239.1", "GL000240.1", "GL000241.1",
            "GL000242.1", "GL000243.1", "GL000244.1", "GL000245.1", "GL000246.1",
            "GL000247.1", "GL000248.1", "GL000249.1", "hs37d5",     "NC_007605"};

        /* 85 compressed mref chromosomes. These are the chromosome sizes + 1. Also, it is unclear,
           why some chromosome sizes differ from the 1K genomes reference, e.g. Chromosome 1 is
           249904550 in there, but significantly smaller here.
           Note that the hardcoded data used to be in MasterMrefProcessor. */
        static const std::vector<ChrSize> chrSizesCompressedMref {
            249250622, 243199374, 198022431, 191154277, 180915261, 171115068,
            159138664, 146364023, 141213432, 135534748, 135006517, 133851896,
            115169879, 107349541, 102531393, 90354754,  81195211,  78077249,
            59128984,  63025521,  48129896,  51304567,  155270561, 59373567,
            106434,    547497,    189790,    191470,    182897,    38915,
            37176,     90086,     169875,    187036,    36149,     40104,
            37499,     81311,     174589,    41002,     4263,      92690,
            159170,    27683,     166567,    186859,    164240,    137719,
            172546,    172295,    172150,    161148,    179199,    161803,
            155398,    186862,    180456,    179694,    211174,    15009,
            128375,    129121,    19914,     43692,     27387,     40653,
            45942,     40532,     34475,     41935,     45868,     39940,
            33825,     41934,     42153,     43524,     43342,     39930,
            36652,     38155,     36423,     39787,     38503,     35477944,
            171824};

        // Used to be -2, but in the mref space 1003 is INVALID,
        // and -2 has the disadvantage that it cannot be represented as an unsigned integer.
        // By moving this to INVALID (1003), we can make CompressedMrefIndex an unsigned integer,
        // and can switch -- for compile-time checks -- the signedness of ChrIndex and
        // CompressedMrefIndex. This gives us a poor-man's type checking, and we can postpone
        // a bigger (more time-consuming) refactoring.
        //
        // Note that NA is only used when mapping from ChrIndex to CompressedMrefIndex, to indicate
        // that chromosome is actually not among the compressed master ref chromosomes.
        static const CompressedMrefIndex NA = 1003;

        // This used to be `indexConverter`.
        static const std::vector<CompressedMrefIndex> indexToCompressedMrefIndex {
            NA, 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17,
            18, 19, 20, 21, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, 22, 23, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
            61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
            80, 81, 82, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 83, 84, NA, NA, NA};

    } /* namespace hg37 */

    const std::string Hg37ChrConverter::assemblyName = "hg37";

    bool Hg37ChrConverter::isValid(ChrIndex index) {
        return index != hg37::INVALID && index != hg37::ZERO && (
            _isAutosome(index) ||
            _isGonosome(index) ||
            _isTechnical(index) ||
            _isVirus(index) ||
            _isExtrachromosomal(index) ||
            _isDecoy(index) ||
            _isUnassigned(index) /* ||  // There are no HLA and ALT contigs in hg37. The ranges are empty.
            _isHLA(index) ||
            _isALT(index) */
        );
    }

    void Hg37ChrConverter::assertValid(ChrIndex index) {
    #ifndef NDEBUG
        if (!isValid(index)) {
            throw_with_trace(std::runtime_error("Invalid chromosome index: " +
                             std::to_string(index)));
        }
    #endif
    }

    bool Hg37ChrConverter::isValid(CompressedMrefIndex index) {
        return index != hg37::NA;
    }

    void Hg37ChrConverter::assertValid(CompressedMrefIndex index) {
    #ifndef NDEBUG
        if (!isValid(index)) {
            throw_with_trace(std::runtime_error("Invalid compressed mref index: " +
                                     std::to_string(index)));
        }
    #endif
    }

    std::vector<ChrIndex> Hg37ChrConverter::_buildCompressedMrefIndexToIndex(
        CompressedMrefIndex nCompressed,
        const std::vector<CompressedMrefIndex> &indexToCompressedMrefIndex) {

        // This is now the only place, where invalid values are assigned, ...
        std::vector<ChrIndex> result (static_cast<unsigned int>(nCompressed), hg37::NA);
        for (ChrIndex globalIndex = 0;
             globalIndex < ChrIndex(indexToCompressedMrefIndex.size());
             ++globalIndex) {

            CompressedMrefIndex compressedMrefIndex =
                indexToCompressedMrefIndex[static_cast<unsigned int>(globalIndex)];

            if (isValid(compressedMrefIndex)) {
                unsigned int cIdx = static_cast<unsigned int>(compressedMrefIndex);
                if (isValid(result[cIdx])) {
                    throw_with_trace(std::runtime_error(
                        "Compressed mref index " + std::to_string(compressedMrefIndex) +
                        " is already assigned to " +
                        std::to_string(result[cIdx]) +
                        " and cannot be assigned to " + std::to_string(globalIndex)));
                }
                result[cIdx] = globalIndex;
            }
        }

        // ... but before we continue, we ensure there are no gaps. There must be an index
        // in the global index space for all compressed mref indices/chromosomes.
        for (auto it = result.cbegin(); it != result.cend(); ++it) {
            assertValid(*it);
        }

        return result;
    }

    Hg37ChrConverter::Hg37ChrConverter(const std::vector<ChrName> &indexToChrName,
                                       const std::vector<ChrName> &compressedMrefIndexToChrName,
                                       const std::vector<ChrSize> &chrSizesCompressedMref,
                                       const std::vector<CompressedMrefIndex> &indexToCompressedMrefIndex) :
                    _indexToChrName {indexToChrName},
                    _compressedMrefIndexToChrName {compressedMrefIndexToChrName},
                    _chrSizesCompressedMref {chrSizesCompressedMref},
                    _indexToCompressedMrefIndex {indexToCompressedMrefIndex},
                    _compressedMrefIndexToIndex {_buildCompressedMrefIndexToIndex(
                        compressedMrefIndexToChrName.size(),
                        indexToCompressedMrefIndex)}{
            if (indexToChrName.size() != indexToCompressedMrefIndex.size())
                throw_with_trace(std::invalid_argument(
                    "indexToChrName and indexToCompressedMrefIndex must have the same size"));
            if (compressedMrefIndexToChrName.size() != chrSizesCompressedMref.size())
                throw_with_trace(std::invalid_argument(
                    "compressedMrefIndexToChrName and chrSizesCompressedMref must have the same size"));
        }

    Hg37ChrConverter::Hg37ChrConverter()
        : Hg37ChrConverter(hg37::indexToChrName,
                           hg37::compressedMrefIndexToChrName,
                           hg37::chrSizesCompressedMref,
                           hg37::indexToCompressedMrefIndex) {}

    ChrIndex Hg37ChrConverter::nChromosomes() const {
        return ChrIndex(_indexToChrName.size());
    }

    CompressedMrefIndex Hg37ChrConverter::nChromosomesCompressedMref() const {
        return CompressedMrefIndex(_compressedMrefIndexToChrName.size());
    }

    /** Map an index position to a chromosome name. */
    ChrName Hg37ChrConverter::indexToChrName(ChrIndex index) const {
//        assertValid(index);
        return _indexToChrName[static_cast<unsigned int>(index)];
    }

    /** chr1-chr22, ... */
    bool Hg37ChrConverter::_isAutosome(ChrIndex index) {
        return hg37::automoseRange.contains(index);
    }
    bool Hg37ChrConverter::isAutosome(ChrIndex index) const {
        return _isAutosome(index);
    }

    /** chrX, chrY */
    bool Hg37ChrConverter::_isGonosome(ChrIndex index) {
        return hg37::gonosomeRange.contains(index);
    }
    bool Hg37ChrConverter::isGonosome(ChrIndex index) const {
        return _isGonosome(index);
    }


    /** phix index. */
    bool Hg37ChrConverter::_isTechnical(ChrIndex index) {
        return hg37::technicalRange.contains(index);
    }
    bool Hg37ChrConverter::isTechnical(ChrIndex index) const {
        return _isTechnical(index);
    }

    /** NC_007605. */
    bool Hg37ChrConverter::_isVirus(ChrIndex index) {
        return hg37::virusRange.contains(index);
    }
    bool Hg37ChrConverter::isVirus(ChrIndex index) const {
        return _isVirus(index);
    }

    /** Mitochondrial chromosome index. */
    bool Hg37ChrConverter::_isExtrachromosomal(ChrIndex index) {
        return hg37::extrachromosomalRange.contains(index);
    }
    bool Hg37ChrConverter::isExtrachromosomal(ChrIndex index) const {
        return _isExtrachromosomal(index);
    }

    /** Decoy sequence index. */
    bool Hg37ChrConverter::_isDecoy(ChrIndex index) {
        return hg37::decoyRange.contains(index);
    }
    bool Hg37ChrConverter::isDecoy(ChrIndex index) const {
        return _isDecoy(index);
    }

    bool Hg37ChrConverter::_isUnassigned(ChrIndex index) {
        return hg37::unassignedRange.contains(index);
    }
    bool Hg37ChrConverter::isUnassigned(ChrIndex index) const {
        return _isUnassigned(index);
    }

    bool Hg37ChrConverter::_isALT(ChrIndex index [[gnu::unused]]) {
        return false;
    }
    bool Hg37ChrConverter::isALT(ChrIndex index [[gnu::unused]]) const {
        return _isALT(index);
    }

    bool Hg37ChrConverter::_isHLA(ChrIndex index [[gnu::unused]]) {
        return false;
    }
    bool Hg37ChrConverter::isHLA(ChrIndex index [[gnu::unused]]) const {
        return _isHLA(index);
    }


    /* Compressed Master Ref chromosomes are 1-22, X, Y, GL* (unassigned), hs37d4 (decoys), and
     * NC_007605 (virus). Excluded are MT and phix. Used to be index <= 1000 (virus). */
    bool Hg37ChrConverter::isCompressedMref(ChrIndex index) const {
//        assertValid(index);
        return isValid(_indexToCompressedMrefIndex.at(static_cast<unsigned int>(index)));
    }

    /** Map an compressed mref index to a chromosome name. */
    ChrName
    Hg37ChrConverter::compressedMrefIndexToChrName(CompressedMrefIndex index) const {
//        assertValid(index);
        return _compressedMrefIndexToChrName.at(static_cast<unsigned int>(index));
    }

    /** Map an index from the global index-space to the compressed mref index-space. */
    CompressedMrefIndex
    Hg37ChrConverter::indexToCompressedMrefIndex(ChrIndex index) const {
//        assertValid(index);
        CompressedMrefIndex result = _indexToCompressedMrefIndex.at(static_cast<unsigned int>(index));
//        assertValid(result);
        return result;
    }

    ChrIndex
    Hg37ChrConverter::compressedMrefIndexToIndex(CompressedMrefIndex index) const {
//        assertValid(index);
        return _compressedMrefIndexToIndex.at(static_cast<unsigned int>(index));
    }

    /** Map compressed mref index to chromosome size. */
    ChrSize
    Hg37ChrConverter::chrSizeCompressedMref(CompressedMrefIndex index) const {
//        assertValid(index);
        return _chrSizesCompressedMref[static_cast<unsigned int>(index)];
    }

    ChrIndex
    Hg37ChrConverter::chrNameToIndex(ChrName chrName) const {
        ChrIndex result;
        try {
            result = parseChrAndReturnIndex(chrName.begin(), chrName.end(), ' ');
        } catch (DomainError &e) {
           throw e << error_info_string("from = " + chrName);
        }
        return result;
    }

    bool
    Hg37ChrConverter::isInBlockedRegion(ChrIndex chrIndex, ChrSize position) const {
//        assertValid(chrIndex);
         // For mate not in range 33140000-33149999 on chromosome 2, do ...
        return !(chrIndex == 2 && (position / 10000 == 3314));
    }

    /* This is parsing code. It takes a position in a character stream, and translates the
       following character(s) into index positions (see ChrConverter::indexToChrName). It is slightly
       modified from the original implementation by Umut Toprak.

       If the first position is a digit, read up to the next stopChar.

         * (\d+)$ -> $1

       If the first position is *not* a digit return indices according to the following rules:

         * h -> 999
         * X -> 40
         * Y -> 41
         * MT -> 1001
         * G?(\d+)\. -> $1
         * N -> 1000
         * p -> 1002

       NOTE: Most of the matches are eager matches, which means the algorithm does not check for
             whether the end iterator or the stopChar is actually reached or whether it follows
             any expected pattern! The actual stopChar is not actually checked in these cases.

       All identifiers not matching any of these rules will throw an exception (domain_error).

       IMPORTANT: The hg37 parser here ignores the stopCharExt, but instead remains with the legacy
                  behavior.
    */
    ChrIndex Hg37ChrConverter::parseChrAndReturnIndex(
            std::string::const_iterator start,
            std::string::const_iterator end,
            char stopChar,
            const std::string &stopCharExt[[gnu::unused]]  // Attribute to remove the warning
            ) const {
        ChrIndex chrIndex {0};
        /* if (start == end) {
            throw_with_trace(DomainError("Chromosome identifier is empty."));
        } else */ if (isdigit(*start)) {
            for (auto chr_cit = start; chr_cit != end && *chr_cit != stopChar; ++chr_cit) {
                chrIndex = chrIndex * 10 + ChrIndex(*chr_cit - '0');
            }
        } else {
            switch (*start) {
                case 'h':
                    chrIndex = hg37::decoyIndex;
                    break;
                case 'X':
                    chrIndex = hg37::xIndex;
                    break;
                case 'G':   // Match GL...... chromosomes.
                    for (auto cit = next(start, 2); *cit != '.'; ++cit) {
                        chrIndex = 10 * chrIndex + ChrIndex(*cit - '0');
                    }
                    break;
                case 'Y':
                    chrIndex = hg37::yIndex;
                    break;
                case 'M':   // Match "MT"
                    ++start;
                    if (start != end && *start == 'T') {
                        chrIndex = hg37::mtIndex;
                    } else {
                        throw_with_trace(
                            DomainError("Chromosome identifier with invalid prefix 'M" +
                                              std::to_string(*start) + "'."));
                    }
                    break;
                case 'N':
                    chrIndex = hg37::decoyIndex;
                    break;
                case 'p':
                    chrIndex = hg37::phixIndex;
                    break;
                default:
                    throw_with_trace(DomainError("Chromosome identifier with invalid prefix '"
                                                       + std::to_string(*start) + "'."));
            }
        }
        return chrIndex;
    }

} /* namespace sophia */
