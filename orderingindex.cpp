#include "stdafx.h"
#include <vector>
#include <list>
#include <queue>
#include <map>
#include <algorithm>
#include <iostream>
#include <deque>
#include <math.h>
#include <time.h>
#include <functional>  
#include <memory>
#include <set>
#include <map>
#include <queue>
#include <string>
#include <memory>
#include <functional>
#include <math.h>
#include <array>
#include <tuple>
using namespace std;
typedef unsigned long long ull;
#pragma region baseinterface
struct IndexRange {
	ull lower;
	ull upper;
	bool contained;
	IndexRange()
	{

	}
	IndexRange(ull a, ull b, bool c =true)
		:lower(a), upper(b), contained(c)
	{

	}
	IndexRange(const IndexRange& o)
	{
		lower = o.lower;
		upper = o.upper;
		contained = o.contained;
	}

	int operator ()(const IndexRange& o)const
	{
		if (lower != o.lower) return lower > o.lower ? 1 : -1;
		if (upper != o.upper) return upper > o.upper ? 1 : -1;

		return 0;
	}
};
/// 用于 std::sort 对结构体按指定成员比较
struct {
	bool operator()(const IndexRange &a, const IndexRange &b) const
	{
		if (a.lower != b.lower) return a.lower > b.lower ? 1 : -1;
		if (a.upper != b.upper) return a.upper > b.upper ? 1 : -1;

		return false;
	}
}CompA;
//也不清楚scala的seq能对应std哪个容器,先用定义后面直接改这个定义即可
using SeqIndexRange = std::list<IndexRange>;
using dblbox = std::array<double, 4>;
using SeqBoxs = std::vector<std::array<double, 4>>;
using SeqTimes = std::vector<std::array<ull, 2>>;
struct XZSFC {
	short DefaultPrecision = 12;
	double LogPointFive = log(0.5);
};
struct ZRange {
	ull min = 0, max = 0, mid = 0, length = 0;

	ZRange(ull umin, ull umax)
	{
		min = umin;
		max = umax;
		mid = (max + min) >> 1;
		length = max - min + 1;
	}


	ZRange(const ZRange& o)
	{
		min = o.min;
		max = o.max;
		mid = o.mid;
		length = o.length;
	}
	bool contains(ull bits)
	{
		return   bits >= min && bits <= max;
	}
	bool contains(const ZRange& r)
	{
		return contains(r.min) && contains(r.max);
	}
	bool overlaps(const ZRange& r) {
		return  contains(r.min) || contains(r.max);
	}
};
using ArrayZRange = std::vector<ZRange>;
class SpaceTimeFillingCurve {
	int FullPrecision = 64;
	int maxRanges = -1;
public:
	virtual ull index(double x, double y, ull t, bool   lenient = false) = 0;
	virtual void invert(ull  z, double& x, double & y, ull&t) = 0;
	virtual  SeqIndexRange ranges(SeqBoxs seqboxs, SeqTimes seqtimes, int precision, int maxRanges) = 0;

	virtual SeqIndexRange ranges(double x1, double x2, double y1, double y2, ull t1, ull t2)
	{
		SeqBoxs d;
		d.resize(1);
		d[0][0] = x1; d[0][1] = y1; d[0][2] = x2; d[0][3] = y1;
		SeqTimes tseq;
		tseq.resize(1);
		tseq[0][0] = (t1);		tseq[0][1] = (t2);
		ranges(d, tseq, FullPrecision, maxRanges);
	}
	virtual SeqIndexRange ranges(double x1, double x2, double y1, double y2, ull t1, ull t2, int precision)
	{
		SeqBoxs d;
		d.resize(1);
		d[0][0] = x1; d[0][1] = y1; d[0][2] = x2; d[0][3] = y1;
		SeqTimes tseq;
		tseq.resize(1);
		tseq[0][0] = (t1);		tseq[0][1] = (t2);
		ranges(d, tseq, precision, maxRanges);
	}
	virtual SeqIndexRange ranges(double x1, double x2, double y1, double y2, ull t1, ull t2, int precision, int maxRanges)
	{
		SeqBoxs d;
		d.resize(1);
		d[0][0] = x1; d[0][1] = y1; d[0][2] = x2; d[0][3] = y1;
		SeqTimes tseq;
		tseq.resize(1);
		tseq[0][0] = (t1);		tseq[0][1] = (t2);
		ranges(d, tseq, precision, maxRanges);
	}

};
class SpaceFillingCurve {
	int FullPrecision = 64;
	int maxRanges = -1;

	virtual ull index(double x, double y, bool   lenient = false) = 0;
	virtual void invert(ull  z, double& x, double & y) = 0;
	virtual  SeqIndexRange ranges(SeqBoxs seqboxs, int precision, int maxRanges) = 0;

	virtual SeqIndexRange ranges(double x1, double x2, double y1, double y2)
	{
		SeqBoxs d;
		d.resize(1);
		d[0][0] = x1; d[0][1] = y1; d[0][2] = x2; d[0][3] = y1;
		ranges(d, FullPrecision, maxRanges);
	}
	virtual SeqIndexRange ranges(double x1, double x2, double y1, double y2, int precision)
	{
		SeqBoxs d;
		d.resize(1);
		d[0][0] = x1; d[0][1] = y1; d[0][2] = x2; d[0][3] = y1;
		ranges(d, precision, maxRanges);
	}
	virtual SeqIndexRange ranges(double x1, double x2, double y1, double y2, int precision, int maxRanges)
	{
		SeqBoxs d;
		d.resize(1);
		d[0][0] = x1; d[0][1] = y1; d[0][2] = x2; d[0][3] = y1;
		ranges(d, precision, maxRanges);
	}

};


struct QueryWindow
{
	double xmin; double xmax; double ymin; double ymax;
	QueryWindow()
	{

	}
	QueryWindow(double dxmin, double dymin, double dxmax, double dymax) :
		xmin(dxmin), xmax(dxmax), ymin(dymin), ymax(dymax)
	{

	}

	QueryWindow(const QueryWindow& o)
	{
		xmin = o.xmin;
		ymin = o.ymin;
		xmax = o.xmax;
		ymax = o.ymax;
	}
};
struct XElement
{
	double xmin; double xmax; double ymin; double ymax; double length;

	double xext = 0;
	double yext = 0;
	XElement() {}
	XElement(double dxmin, double dymin, double dxmax, double dymax, double dlength)
		:xmin(dxmin), xmax(dxmax), ymin(dymin), ymax(dymax), length(dlength)
	{
		xext = xmax + length;
		yext = ymax + length;
	}
	XElement(const XElement& o)
	{
		xmin = o.xmin;
		xmax = o.xmax;
		ymin = o.ymin;
		ymax = o.ymax;
		length = o.length;
		xext = o.xext;
		yext = o.yext;
	}
	XElement clone()
	{
		return XElement(xmin, ymin, xmax, ymax, length);
	}
	bool isContained(QueryWindow window)
	{
		return window.xmin <= xmin && window.ymin <= ymin && window.xmax >= xext && window.ymax >= yext;
	}
	bool overlaps(QueryWindow window)
	{
		return 	window.xmax >= xmin && window.ymax >= ymin && window.xmin <= xext && window.ymin <= yext;
	}
	std::vector<XElement*> children() {
		double xCenter = (xmin + xmax) / 2.0;
		double yCenter = (ymin + ymax) / 2.0;
		double len = length / 2.0;
		std::vector<XElement*> ret;
		XElement* c0 = new XElement(xmin, ymin, xmax, ymax, length);
		c0->xmax = xCenter;
		c0->ymax = yCenter;
		c0->length = len;
		XElement* c1 = new XElement(xmin, ymin, xmax, ymax, length);
		c1->xmin = xCenter; c1->ymax = yCenter; c1->length = len;
		XElement* c2 = new XElement(xmin, ymin, xmax, ymax, length);
		c2->xmax = xCenter; c2->ymin = yCenter; c2->length = len;
		XElement* c3 = new XElement(xmin, ymin, xmax, ymax, length);
		c3->xmin = xCenter; c3->ymin = yCenter; c3->length = len;
		ret.push_back(c0);
		ret.push_back(c1);
		ret.push_back(c2);
		ret.push_back(c3);
		return ret;
	}
};
struct BitNormalizedDimension {
	double min, max;
	int maxIndex;
	double normalizer;
	double denormalizer;
	BitNormalizedDimension()
	{

	}
	BitNormalizedDimension(const BitNormalizedDimension & o)
	{
		min = o.min;
		max = o.max;
		maxIndex = o.maxIndex;
		normalizer = o.denormalizer;
		denormalizer = o.denormalizer;
	}
	BitNormalizedDimension(double dmin, double dmax, int precision)
	{
		if (!(precision > 0 && precision < 32))
			return;
		min = dmin;
		max = dmax;
		ull bins = pow(2, precision);// .toLong1 << precision;
		normalizer = bins / (max - min);
		denormalizer = (max - min) / bins;
		maxIndex = (bins - 1);
	}


	//case class NormalizedLat(precision : Int) extends BitNormalizedDimension(-90d, 90d, precision)
	//case class NormalizedLon(precision : Int) extends BitNormalizedDimension(-180d, 180d, precision)
	//case class NormalizedTime(precision : Int, override val max : Double) extends BitNormalizedDimension(0d, max, precision)


	int normalize(double x) {
		if (x >= max)
			return maxIndex;
		else
			return  floor((x - min) * normalizer);
	}
	double denormalize(int x) {
		if (x >= maxIndex)
		{
			return min + (maxIndex + 0.5) * denormalizer;
		}
		else
		{
			return min + (x + 0.5) * denormalizer;
		}
	}
};
#pragma endregion

#pragma region ZNOrdering
struct ZPrefix
{
	ull prefix;
	int precision;
	ZPrefix()
	{
	}
	ZPrefix(ull llprefix, int iprecision)
	{
		prefix = llprefix;
		precision = iprecision;
	}
	ZPrefix(ZPrefix&o)
	{
		prefix = o.prefix;
		precision = o.precision;
	}
};
/**
  * N-dimensional z-curve base class
  */
struct ZN {

	// number of bits used to store each dimension
	int BitsPerDimension;

	// number of dimensions
	int Dimensions;

	// max value for this z object - can be used to mask another long using &
	ull MaxMask;

	// total bits used - usually bitsPerDim * dims
	int TotalBits;

	// number of quadrants in our quad/oct tree - has to be lazy to instantiate correctly
	double Quadrants;
	int DefaultRecurse = 7;
	ZN(int Dimensions)
	{
		Quadrants = pow(2, Dimensions);
	}
	/**
	  * Insert (Dimensions - 1) zeros between each bit to create a zvalue from a single dimension.
	  * Only the first BitsPerDimension can be considered.
	  *
	  * @param value value to split
	  * @return
	  */
	virtual ull split(ull value) = 0;

	/**
	  * Combine every (Dimensions - 1) bits to re-create a single dimension. Opposite of split.
	  *
	  * @param z value to combine
	  * @return
	  */
	virtual int combine(ull z) = 0;

	/**
	  * Is the value contained in the range. Considers user-space.
	  *
	  * @param range range
	  * @param value value to be tested
	  * @return
	  */
	virtual bool contains(ZRange &range, ull value) = 0;

	/**
	  * Is the value contained in the range. Considers user-space.
	  *
	  * @param range range
	  * @param value value to be tested
	  * @return
	  */
	bool contains(ZRange& range, ZRange value)
	{
		return contains(range, value.min) && contains(range, value.max);
	}

	/**
	  * Does the value overlap with the range. Considers user-space.
	  *
	  * @param range range
	  * @param value value to be tested
	  * @return
	  */
	virtual bool overlaps(ZRange& range, ZRange& value) = 0;

	/**
	  * Returns (litmax, bigmin) for the given range and point
	  *
	  * @param p point
	  * @param rmin minimum value
	  * @param rmax maximum value
	  * @return (litmax, bigmin)
	  */
	void zdivide(ull p, ull  rmin, ull rmax, ull& minz, ull&maxz)
	{
		//: (Long, Long) = ZN.zdiv(load, Dimensions)(p, rmin, rmax)
	}
	// virtual  SeqIndexRange ranges(SeqBoxs seqboxs, int precision, int maxRanges) = 0;
	virtual SeqIndexRange zranges(ArrayZRange & zbounds)
	{
		return  zranges(zbounds, 64, -1, DefaultRecurse);
	}
	virtual SeqIndexRange zranges(ArrayZRange & zbounds, int precision)
	{
		return  zranges(zbounds, precision, -1, DefaultRecurse);
	}
	virtual SeqIndexRange zranges(ArrayZRange & zbounds, int precision, int maxRanges)
	{
		return  zranges(zbounds, precision, maxRanges, DefaultRecurse);
	}
	/**
* Calculates ranges in index space that match any of the input bounds. Uses breadth-first searching to
* allow a limit on the number of ranges returned.
*
* To improve performance, the following decisions have been made:
*   uses loops instead of foreach/maps
*   uses java queues instead of scala queues
*   allocates initial sequences of decent size
*   sorts once at the end before merging
*
* @param zbounds search space
* @param precision precision to consider, in bits (max 64)
* @param maxRanges loose cap on the number of ranges to return. A higher number of ranges will have less
*                  false positives, but require more processing.
* @param maxRecurse max levels of recursion to apply before stopping
* @return ranges covering the search space
*/
	struct   ZPrefix
	{
		ull prefix;
		int precision;
	};
	virtual SeqIndexRange zranges(ArrayZRange & zbounds, int precision = 64, int maxRanges = -1, int maxRecurse = 7)
	{

		// stores our results - initial size of 100 in general saves us some re-allocation
		//val ranges = new java.util.ArrayList[IndexRange](100)
		std::list<IndexRange> ranges;
		ranges.resize(100);
		// values remaining to process - initial size of 100 in general saves us some re-allocation
		//val remaining = new java.util.ArrayDeque[(Long, Long)](100)
		std::deque<std::pair<ull, ull>> remaining;
		remaining.resize(100);
		// calculate the common prefix in the z-values - we start processing with the first diff
		ZPrefix tmpZPrefix = longestCommonPrefix(zbounds);
		ull commonPrefix = tmpZPrefix.prefix;
		int	commonBits = tmpZPrefix.precision;
		int offset = 64 - commonBits;

		// checks if a range is contained in the search space

		auto isContained = [&zbounds, this](ZRange& range)
		{
			int i = 0;
			while (i < zbounds.size())
			{
				if (contains(zbounds[i], range))
					return true;
				i++;
			}
			return false;

		};
		auto isOverlapped = [&zbounds, this](ZRange& range)
		{
			int i = 0;
			while (i < zbounds.size())
			{
				if (overlaps(zbounds[i], range))
					return true;
				i++;
			}
			return false;

		};

		auto  checkValue = [&offset, &isContained, &isOverlapped, &precision, &ranges, &remaining](ull prefix, ull  quadrant) ->unsigned int {
			ull min = prefix | (quadrant << offset);
			ull max = min | (1L << offset) - 1;
			ZRange quadrantRange = ZRange(min, max);

			if (isContained(quadrantRange) || offset < 64 - precision) {
				// whole range matches, happy day
				ranges.push_back(IndexRange(min, max, true));
			}
			else if (isOverlapped(quadrantRange)) {
				// some portion of this range is excluded
				// queue up each sub-range for processing
				remaining.push_back(std::make_pair(min, max));
			}
		};
		std::pair<ull, ull> LevelTerminator = std::make_pair(-1, -1);
		auto  bottomOut = [this, &remaining, &ranges, LevelTerminator]() ->unsigned int {
			do {
				std::pair<ull, ull> minMax = remaining.front();
				auto min2 = std::get<0>(LevelTerminator);
				auto max2 = std::get<1>(LevelTerminator);
				if (minMax.first == min2 && minMax.second == max2)
				{
					ranges.emplace_back(minMax.first, minMax.second, false);
				}
				remaining.pop_front();
			} while (!remaining.empty());
		};


		// initial level - we just check the single quadrant
		checkValue(commonPrefix, 0);
		remaining.push_back(LevelTerminator);
		offset -= Dimensions;

		// level of recursion
		int level = 0;

		int rangeStop = maxRanges == -1 ? INT_MAX : maxRanges;
		int recurseStop = (maxRecurse == 7) ? (DefaultRecurse) : maxRecurse;

		do {
			auto next = remaining.front();
			remaining.pop_front();
			if (next.first == LevelTerminator.first && next.second == LevelTerminator.second)
			{
				if (!remaining.empty()) {
					level += 1;
					offset -= Dimensions;
					if (level >= recurseStop || offset < 0) {
						bottomOut();
					}
					else {
						remaining.push_back(LevelTerminator);
					}
				}
			}
			else {
				auto prefix = next.first;
				ull quadrant = 0L;
				while (quadrant < Quadrants) {
					checkValue(prefix, quadrant);
					quadrant += 1;
				};
				// subtract one from remaining.size to account for the LevelTerminator
				if (ranges.size() + remaining.size() - 1 >= rangeStop) {
					bottomOut();
				}
			}
		} while (!remaining.empty());

		// we've got all our ranges - now reduce them down by merging overlapping values
		//std::sort(ranges.begin(), ranges.end());
		ranges.sort( CompA);

		auto current = ranges.front(); // note: should always be at least one range
		SeqIndexRange result;
		int i = 1;
		while (i < ranges.size()) {
			auto range = ranges.front();
			if (range.lower <= current.upper + 1) {
				// merge the two ranges
				current = IndexRange(current.lower, max(current.upper, range.upper), current.contained && range.contained);
			}
			else {
				// append the last range and set the current range for future merging
				result.push_back(current);
				current = range;
			}
			i += 1;
		}
		// append the last range - there will always be one left that wasn't added
		result.push_back(current);

		return		  result;
	}

	/**
	 * Cuts Z-Range in two and trims based on user space, can be used to perform augmented binary search
	 *
	 * @param xd: division point
	 * @param inRange: is xd in query range
	 */
	static SeqIndexRange g_empty;
	SeqIndexRange cut(ZRange &r, ull xd, bool inRange) {
		SeqIndexRange ret;
		if (r.min == r.max) {
			return g_empty;
		}
		else if (inRange) {
			if (xd == r.min) { // degenerate case, two nodes min has already been counted
				ret.emplace_back(r.max, r.max);
			}
			else if (xd == r.max) { // degenerate case, two nodes max has already been counted
				ret.emplace_back(r.min, r.min);
			}
			else {
				ret.emplace_back(r.min, xd - 1);
				ret.emplace_back(xd + 1, r.max);
			}
		}
		else {
			ull litmax0, bigmin = 0;
			std::tuple<ull, ull> tup = zdiv(Dimensions, xd, r.min, r.max);

			ret.emplace_back(r.min, std::get<0>(tup));
			ret.emplace_back(std::get<1>(tup), r.max);
		}
	}

	/**
	  * Calculates the longest common binary prefix between two z longs
	  *
	  * @return (common prefix, number of bits in common)
	  */
	  //这里不是指针是不确定参数大小
	ZPrefix longestCommonPrefix(ArrayZRange & values) {
		if (values.size() < 1)
			return ZPrefix();
		std::vector<ull> llu;
		for (auto& var : values)
		{
			llu.push_back(var.min);
			llu.push_back(var.max);
		}

		int bitShift = TotalBits - Dimensions;
		ull head = llu[0] >> bitShift;
		int tail = llu.size() - 1;
		for (int i = llu.size() - 1; i > 0; i--)
		{
			ull v = llu[i] >> bitShift;
			if ((v == head) && bitShift > -1)
			{
				bitShift -= Dimensions;
				head = llu[0] >> bitShift;
			}
		}

		bitShift += Dimensions; // increment back to the last valid value
		ZPrefix ret;
		ret.prefix = (llu[0] & (((ull)~((ull)0)) << bitShift));
		ret.precision = 64 - bitShift;
		return ret;
	}
	//	case class ZPrefix(prefix : Long, precision : Int) // precision in bits
		/** Loads either 1000... or 0111... into starting at given bit index of a given dimension */
	ull load(ull target, ull p, int bits, int dim)
	{
		auto mask = ~(split(MaxMask >> (BitsPerDimension - bits)) << dim);
		auto wiped = target & mask;
		return wiped | (split(p) << dim);
	}





	/**
	  * Implements the the algorithm defined in: Tropf paper to find:
	  * LITMAX: maximum z-index in query range smaller than current point, xd
	  * BIGMIN: minimum z-index in query range greater than current point, xd
	  *
	  * @param load: function that knows how to load bits into appropraite dimension of a z-index
	  * @param xd: z-index that is outside of the query range
	  * @param rmin: minimum z-index of the query range, inclusive
	  * @param rmax: maximum z-index of the query range, inclusive
	  * @return (LITMAX, BIGMIN)
	  */
private:
	//std::tuple<ull, ull> LevelTerminator = std::make_tuple(-1,-1);

	std::tuple<ull, ull> zdiv(int dims, ull xd, ull  rmin, ull rmax)
	{
		//require(rmin < rmax, "min ($rmin) must be less than max $(rmax)")
		ull zmin = rmin;
		ull zmax = rmax;
		ull bigmin = 0L;
		ull litmax = 0L;
		//int dims = 2;
		auto bit = [&](ull x, int idx) {	return (int)((x & (1ULL << idx)) >> idx);	};
		auto over = [&](ull bits) {return  1L << (bits - 1); };
		auto under = [&](ull bits) { return  (1L << (bits - 1)) - 1; };


		auto i = 64;
		while (i > 0) {
			i -= 1;

			auto bits = i / dims + 1;
			auto dim = i % dims;
			auto a = bit(xd, i);
			auto b = bit(zmin, i);
			auto c = bit(zmax, i);
			if (a == 0 && b == 0 && c == 1)
			{
				zmax = load(zmax, under(bits), bits, dim);
				bigmin = load(zmin, over(bits), bits, dim);
			}
			else if (a == 0 && b == 1 && c == 1)
			{
				bigmin = zmin;
				return std::make_tuple(litmax, bigmin);
			}

			else if (a == 1 && b == 0 && c == 0)
			{
				litmax = zmax;
				return std::make_tuple(litmax, bigmin);
			}
			else if (a == 1 && b == 0 && c == 1)
			{
				litmax = load(zmax, under(bits), bits, dim);
				zmin = load(zmin, over(bits), bits, dim);
			}
			else
			{
				continue;
			}
		};
		return std::make_tuple(litmax, bigmin);
	}
};



#pragma endregion


#pragma region ZOrdering  = Z2SFC.scala


struct Z2:public ZN
{
	int Dimensions = 2;
	int BitsPerDimension = 31;
	int TotalBits = 62;
	ull MaxMask = 0x7fffffffL;
	ull z = 0;
	int d0 = 0, d1 = 0;
	int DefaultRecurse = 7;
	Z2(int dim = 2): ZN(2)
	{

	}
	Z2(ull Z) :ZN(2),z(Z)
	{
		d0 = dim(0);
		d1 = dim(1);
	}
	Z2(int x, int y):ZN(2)
	{
		z = split(x) | (split(y) << 1);
		d0 = dim(0);
		d1 = dim(1);
	}
	//def <(other: Z2) = z < other.z
	//def >(other: Z2) = z > other.z
	//def <= (other: Z2) = z <= other.z
	//def >= (other: Z2) = z >= other.z
	//def + (offset: Long) = new Z2(z + offset)
	//def - (offset: Long) = new Z2(z - offset)
	//def == (other: Z2) = other.z == z

	virtual int combine(ull z)
	{
		ull x = z & 0x5555555555555555L;
		x = (x ^ (x >> 1)) & 0x3333333333333333L;
		x = (x ^ (x >> 2)) & 0x0f0f0f0f0f0f0f0fL;
		x = (x ^ (x >> 4)) & 0x00ff00ff00ff00ffL;
		x = (x ^ (x >> 8)) & 0x0000ffff0000ffffL;
		x = (x ^ (x >> 16)) & 0x00000000ffffffffL;
		return x;
	}
	ull split(ull value)
	{
		ull x = value & MaxMask;
		x = (x ^ (x << 32)) & 0x00000000ffffffffL;
		x = (x ^ (x << 16)) & 0x0000ffff0000ffffL;
		x = (x ^ (x << 8)) & 0x00ff00ff00ff00ffL; // 11111111000000001111111100000000..
		x = (x ^ (x << 4)) & 0x0f0f0f0f0f0f0f0fL; // 1111000011110000
		x = (x ^ (x << 2)) & 0x3333333333333333L; // 11001100..
		x = (x ^ (x << 1)) & 0x5555555555555555L;// 1010...
		return x;
	}
	void decode(int &x, int &y)
	{
		x = combine(z);
		y = combine(z >> 1);
	}
	int dim(int i) {
		return   Z2::combine(z >> i);
	}
	Z2	  mid(Z2& p)
	{
		int x, y;
		decode(x, y);
		int px, py;
		p.decode(px, py);
		return   Z2((x + px) >> 1, (y + py) >> 1); // overflow safe mean
	}
	//包含
	virtual bool  contains(ZRange & range, ull value)
	{
		Z2 z2tmp(value);
		int x, y;
		z2tmp.decode(x, y);
		return x >= Z2(range.min).d0 && x <= Z2(range.max).d0 && y >= Z2(range.min).d1 && y <= Z2(range.max).d1;
	}
	virtual bool overlaps(int a1, int a2, int b1, int b2)
	{
		return  max(a1, b1) <= min(a2, b2);
	}
	virtual bool  overlaps(ZRange& range, ZRange& value)
	{
		return	overlaps(Z2(range.min).d0, Z2(range.max).d0, Z2(value.min).d0, Z2(value.max).d0) &&
			overlaps(Z2(range.min).d1, Z2(range.max).d1, Z2(value.min).d1, Z2(value.max).d1);
	}


};

struct Z2SpaceFillingCurve : public SpaceFillingCurve
{
	short DefaultPrecision = 31;
	double LogPointFive = log(0.5);

	BitNormalizedDimension lon;
	BitNormalizedDimension lat;
	Z2SpaceFillingCurve()
	{

	}
	Z2SpaceFillingCurve(int precision)
	{
		lon = BitNormalizedDimension(-180., 180., precision);
		lat = BitNormalizedDimension(-180.0, 180.0, precision);
	}
	//正算
	ull  index(double x, double y, bool lenient = false)
	{
		if (lenient)
		{
			return lenientIndex(x, y);
		}
		return Z2(lon.normalize(x), lat.normalize(y)).z;
	}
	//正算
	ull lenientIndex(double x, double y)
	{
		double bx = x;
		if (x < lon.min) { x = lon.min; }
		else if (x > lon.max) { x = lon.max; }

		double by = y;
		if (y < lat.min) { y = lat.min; }
		else if (y > lat.max) { y = lat.max; }

		return Z2(lon.normalize(bx), lat.normalize(by)).z;
	}
	//反算
	void invert(unsigned long long z, double& x, double & y)
	{
		int xr, yc;
		Z2(z).decode(xr, yc);
		x = lon.denormalize(xr);
		y = lat.denormalize(yc);
	}

	virtual  SeqIndexRange ranges(SeqBoxs seqboxs, SeqTimes seqtimes, int precision, int maxRanges)
	{
		ArrayZRange az;
		for (int i = 0; i < seqboxs.size(); i++)
		{
			ull min = index(seqboxs[i][0], seqboxs[i][1]);
			ull max = index(seqboxs[i][2], seqboxs[i][3]);
			az.emplace_back(min, max);
		}
		//z2先写的,这里不改了, scala代码写的真
		Z2 zn(2);
		return zn.zranges(az, precision, maxRanges, 7);
	}
};
#pragma endregion

#pragma region XZOrdering =XZ2SFC.scala
struct XZ2SpaceFillingCurve
{
	double xLo = 0;
	double xHi = 0;
	double yLo = 0;
	double yHi = 0;

	double xSize = 0;
	double  ySize = 0;
	short DefaultPrecision = 12;
	double	 LogPointFive = log(0.5);
	short g = 0;
	XZ2SpaceFillingCurve(short sg, double xmin, double ymin, double xmax, double ymax)
	{
		g = sg;
		xLo = xmin;
		xHi = xmax;
		yLo = ymin;
		yHi = ymax;
		xSize = xHi - xLo;
		ySize = yHi - yLo;
	}
	bool predicate(double min, double max, double w2)
	{
		return max <= (floor(min / w2) * w2) + (2 * w2);
	};
	virtual ull index(double xmin, double ymin, double xmax, double ymax, bool lenient = false)
	{
		// normalize inputs to [0,1]
		double nxmin, nymin, nxmax, nymax;
		normalize(xmin, ymin, xmax, ymax, lenient, nxmin, nymin, nxmax, nymax);

		// calculate the length of the sequence code (section 4.1 of XZ-Ordering paper)

		double maxDim = max(nxmax - nxmin, nymax - nymin);

		// l1 (el-one) is a bit confusing to read, but corresponds with the paper's definitions
		int l1 = floor(log(maxDim) / LogPointFive);

		// the length will either be (l1) or (l1 + 1)
		ull length = 0;
		if (l1 >= g)
		{
			length = g;
		}
		else
		{
			double w2 = pow(0.5, l1 + 1);
			if (predicate(nxmin, nxmax, w2) && predicate(nymin, nymax, w2))
				length = l1 + 1;
			else
				length = l1;
		}

		return sequenceCode(nxmin, nymin, length);
	}
	virtual void invert(ull  z, double& x, double & y, ull&t)
	{
	}
	ull  sequenceCode(double x, double y, int length) {
		double xmin = 0.0;
		double ymin = 0.0;
		double xmax = 1.0;
		double ymax = 1.0;

		ull cs = 0;

		int i = 0;
		while (i < length) {
			double xCenter = (xmin + xmax) / 2.0;
			double yCenter = (ymin + ymax) / 2.0;
			bool t1 = x < xCenter;
			bool t2 = y < yCenter;

			if (t1 == true && t2 == true)
			{
				cs += 1L; xmax = xCenter; ymax = yCenter;
			}
			else if (t1 == false && t2 == true)
			{
				cs += 1L + 1L * ((ull)pow(4, g - i) - 1L) / 3L; xmin = xCenter; ymax = yCenter;
			}
			else if (t1 == true && t2 == false) {
				cs += 1L + 2L * ((ull)pow(4, g - i) - 1L) / 3L; xmax = xCenter; ymin = yCenter;
			}
			else  if (t1 == false && t2 == false)
			{
				cs += 1L + 3L * ((ull)pow(4, g - i) - 1L) / 3L; xmin = xCenter; ymin = yCenter;
			}
			i += 1;
		}

		return cs;
	}

	void sequenceInterval(double x, double y, short length, bool partial, ull& omin, ull&omax)
	{
		ull min = sequenceCode(x, y, length);
		ull max = min;
		if (partial)
		{
			max = min;
		}
		else
		{
			max = min + ((ull)pow(4, g - length + 1) - 1L) / 3L;
		}
		omin = min;
		omax = max;
	}
	void  normalize(double xmin, double ymin, double xmax, double ymax, bool lenient,
		double &dxmin, double&dymin, double &dxmax, double& dymax)
	{
		dxmin = (xmin - xLo) / xSize;
		dymin = (ymin - yLo) / ySize;
		dxmax = (xmax - xLo) / xSize;
		dymax = (ymax - yLo) / ySize;

	}

	// the initial level of quads
	std::vector<XElement*> LevelOneElements = XElement(0.0, 0.0, 1.0, 1.0, 1.0).children();
	XElement* LevelTerminator = new XElement(-1.0, -1.0, -1.0, -1.0, 0);
	std::map<short, XZ2SpaceFillingCurve*> cache;
	XZ2SpaceFillingCurve(short g)
	{
		auto sfc = cache.find(g);
		if (sfc == cache.end())
		{
			XZ2SpaceFillingCurve *sfct = new XZ2SpaceFillingCurve(g, -180.0, -90, 180.0, 90);
			cache[g] = sfct;
		}
	}

	virtual	SeqIndexRange ranges(vector<QueryWindow>& query, int rangeStop = INT_MAX)
	{
		// stores our results - initial size of 100 in general saves us some re-allocation
		vector<IndexRange> ranges;
		//ranges;
		// values remaining to process - initial size of 100 in general saves us some re-allocation
		deque<XElement*> remaining;
		//remaining.resize(100);

		auto isContained = [&query](XElement& quad)
		{
			int i = 0;
			while (i < query.size())
			{
				if (quad.isContained(query[i]))
				{
					return true;
				}
				i += 1;
			};
			return false;
		};


		// checks if a quad overlaps the search space
		auto isOverlapped = [&query](XElement & quad)
		{
			int i = 0;
			while (i < query.size())
			{
				if (quad.overlaps(query[i]))
				{
					return true;
				}
				i += 1;
			}
			return false;
		};

		XZ2SpaceFillingCurve* ths = this;
		auto checkValue = [&ranges, &remaining, isContained, isOverlapped, ths](XElement * quad, short level)
		{
			ull  min, max;
			if (isContained(*quad))
			{
				ths->sequenceInterval(quad->xmin, quad->ymin, level, false, min, max);
				ranges.emplace_back(min, max, true);
			}
			else if (isOverlapped(*quad))
			{
				ths->sequenceInterval(quad->xmin, quad->ymin, level, true, min, max);
				ranges.emplace_back(min, max, false);
				auto cds = quad->children();
				for (int i = 0; i < cds.size(); i++)
					remaining.push_back(cds[i]);

			}
		};


		for (int i = 0; i < LevelOneElements.size(); i++)
			remaining.push_back(LevelOneElements[i]);
		remaining.push_back(LevelTerminator);
		// level of recursion
		short level = 1;

		while (level < g && !remaining.empty() && ranges.size() < rangeStop)
		{
			auto next = remaining.front();
			remaining.pop_front();
			if (next == (LevelTerminator))
			{
				if (!remaining.empty())
				{
					level = (level + 1);
					remaining.push_back(LevelTerminator);
				}
			}
			else {
				checkValue(next, level);
			}
		};

		while (!remaining.empty()) {
			auto quad = remaining.front();
			remaining.pop_front();
			if (quad == LevelTerminator)
			{
				level = (level + 1);
			}
			else {
				ull min, max;
				sequenceInterval(quad->xmin, quad->ymin, level, false, min, max);
				ranges.emplace_back(min, max, false);
			}
		}

		std::sort(ranges.begin(), ranges.end(), CompA);


		auto current = ranges[0];
		list<IndexRange> result;
		int i = 1;
		while (i < ranges.size()) {
			auto range = ranges[i];
			if (range.lower <= current.upper + 1)
			{
				current = IndexRange(current.lower, max(current.upper, range.upper), current.contained && range.contained);
			}
			else {
				// append the last range and set the current range for future merging
				result.push_back(current);
				current = range;
			}
			i += 1;
		}
		result.push_back(current);
		return result;
	}
};
#pragma endregion

#pragma region Z3Ordering =Z3SFC.scala

enum class TimePeriod
{
	eDay,
	eWeek,
	eMonth,
	eYear,
};
struct BinnedTime
{

	//case TimePeriod.Day = > ChronoUnit.DAYS.getDuration.toMillis
	//case TimePeriod.Week = > ChronoUnit.WEEKS.getDuration.toMillis / 1000L
	//case TimePeriod.Month = > (ChronoUnit.DAYS.getDuration.toMillis / 1000L) * 31L
	// based on 365 days + 1 leap day, with a fudge factor of 10 minutes to account for leap seconds added each year
	//case TimePeriod.Year = > (ChronoUnit.DAYS.getDuration.toMinutes * 366L) + 10L
	//86400000
	//604800
	//2678400
	//527050
	static ull maxOffset(TimePeriod period)
	{
		switch (period)
		{
		case TimePeriod::eDay:
			return 86400000;

		case TimePeriod::eWeek:
			return 604800;

		case TimePeriod::eMonth:
			return	2678400;

		case TimePeriod::eYear:
			return 527050;
		default:
			break;
		}

	}
};


struct  Z3 :public ZN{
	int Dimensions = 2;
	int BitsPerDimension = 31;
	int TotalBits = 62;
	ull MaxMask = 0x1fffffL;
	ull z = 0;
	int d0 = 0, d1 = 0, d2 = 0;
	int DefaultRecurse = 7;
	Z3(int dim = 3):ZN(dim)
	{
		Dimensions = dim;
	}
	Z3(ull z):ZN(3)
	{
		d0 = combine(z);
		d1 = combine(z >> 1);
		d2 = combine(z >> 2);
	}
	Z3(int x, int y, int t) :ZN(3)
	{
		ull z = (split(x) | split(y) << 1 | split(t) << 2);
		d0 = combine(z);
		d1 = combine(z >> 1);
		d2 = combine(z >> 2);
	}
	//def < (other: Z3) = z < other.z
	//def >(other: Z3) = z > other.z
	//def >= (other: Z3) = z >= other.z
	//def <= (other: Z3) = z <= other.z
	//def + (offset: Long) = new Z3(z + offset)
	//def - (offset: Long) = new Z3(z - offset)
	//def == (other: Z3) = other.z == z
	//def decode : (Int, Int, Int) = (d0, d1, d2)

	int dim(int i)
	{
		if (i == 0)	return d0;
		else if (i == 1)return  d1;
		else if (i == 2)return   d2;
		else {
			return 0;//throw new IllegalArgumentException(s"Invalid dimension $i - valid dimensions are 0,1,2")
		}
	}

	bool inRange(Z3 &rmin, Z3 &rmax)
	{
		int x = d0; int  y = d1; int z = d2;
		return x >= rmin.d0 &&
			x <= rmax.d0 &&
			y >= rmin.d1 &&
			y <= rmax.d1 &&
			z >= rmin.d2 &&
			z <= rmax.d2;
	}

	Z3 mid(Z3 &p)
	{
		int x = d0; int  y = d1; int z = d2;
		int px = p.d0; int py = p.d1; int pz = p.d2;
		return  Z3((x + px) >> 1, (y + py) >> 1, (z + pz) >> 1);
	}
	//toString 方法未实现
	//def bitsToString = f"(${z.toBinaryString.toLong}%016d)(${d0.toBinaryString.toLong}%08d,${d1.toBinaryString.toLong}%08d,${d2.toBinaryString.toLong}%08d)"
	//override def toString = f"$z $decode"
	//def unapply(z : Z3) : Option[(Int, Int, Int)] = Some(z.decode)
	/** insert 00 between every bit in value. Only first 21 bits can be considered. */
	virtual ull split(ull  value)
	{
		auto x = value & MaxMask;
		x = (x | x << 32) & 0x1f00000000ffffL;
		x = (x | x << 16) & 0x1f0000ff0000ffL;
		x = (x | x << 8) & 0x100f00f00f00f00fL;
		x = (x | x << 4) & 0x10c30c30c30c30c3L;
		return (x | x << 2) & 0x1249249249249249L;
	}

	/** combine every third bit to form a value. Maximum value is 21 bits. */
	virtual int combine(ull  value)
	{
		auto x = z & 0x1249249249249249L;
		x = (x ^ (x >> 2)) & 0x10c30c30c30c30c3L;
		x = (x ^ (x >> 4)) & 0x100f00f00f00f00fL;
		x = (x ^ (x >> 8)) & 0x1f0000ff0000ffL;
		x = (x ^ (x >> 16)) & 0x1f00000000ffffL;
		x = (x ^ (x >> 32)) & MaxMask;
		return (int)x;// x.toInt
	}

	virtual bool contains(ZRange& range, ull value)
	{
		Z3 ztmp(value);
		int x = ztmp.d0, y = ztmp.d1, z = ztmp.d2;
		return x >= Z3(range.min).d0 && x <= Z3(range.max).d0 &&
			y >= Z3(range.min).d1 && y <= Z3(range.max).d1 &&
			z >= Z3(range.min).d2 && z <= Z3(range.max).d2;
	}
	virtual bool overlaps(int a1, int a2, int b1, int b2)
	{
		return  max(a1, b1) <= min(a2, b2);
	}
	virtual bool overlaps(ZRange &range, ZRange& value)
	{
		return overlaps(Z3(range.min).d0, Z3(range.max).d0, Z3(value.min).d0, Z3(value.max).d0) &&
			overlaps(Z3(range.min).d1, Z3(range.max).d1, Z3(value.min).d1, Z3(value.max).d1) &&
			overlaps(Z3(range.min).d2, Z3(range.max).d2, Z3(value.min).d2, Z3(value.max).d2);
	}
};


/**
  * Z3 space filling curve
  *
  * @param period time period used to bin results
  * @param precision bits used per dimension - note all precisions must sum to less than 64
  */
class Z3SpaceTimeFillingCurve :public  SpaceTimeFillingCurve {
	BitNormalizedDimension lat;// (-90, 90, precision);
	BitNormalizedDimension lon;// (-180, 180, precision);

	BitNormalizedDimension time;// (0, maxtime, precision);
public:
	Z3SpaceTimeFillingCurve(TimePeriod period, int precision = 21)
	{
		if (precision <= 0 || precision >= 22)
		{
			printf("error Z3SpaceTimeFillingCurve Precision (bits) per dimension must be in [1,21]");
			return;
		}
		lat = BitNormalizedDimension(-90, 90, precision);
		lon = BitNormalizedDimension(-180, 180, precision);
		double maxtime = BinnedTime::maxOffset(period);
		lon = BitNormalizedDimension(0, maxtime, precision);
		SeqTimes wholePeriod;
		wholePeriod.resize(1);
		wholePeriod[0][0] = time.min;
		wholePeriod[0][1] = time.max;

	}
	virtual ull index(double x, double y, ull t, bool   lenient = false)
	{

		if (x >= lon.min && x <= lon.max && y >= lat.min && y <= lat.max && t >= time.min && t <= time.max)
		{
			return Z3(lon.normalize(x), lat.normalize(y), time.normalize(t)).z;
		}
		else
		{
			return 0;
		}

		//catch {
		//case _: IllegalArgumentException if lenient = > lenientIndex(x, y, t)

	}
	ull lenientIndex(double x, double y, ull  t)
	{
		double bx = 0; if (x < lon.min) { bx = lon.min; }
		else if (x > lon.max) { bx = lon.max; }
		else { bx = x; }
		double by = 0; if (y < lat.min) { by = lat.min; }
		else if (y > lat.max) { by = lat.max; }
		else { by = y; }
		double bt = 0; if (t < time.min) { bt = time.min; }
		else if (t > time.max) { bt = time.max; }
		else { bt = t; }
		return Z3(lon.normalize(bx), lat.normalize(by), time.normalize(bt)).z;
	}
	virtual void invert(ull  z, double& x, double & y, ull&t)
	{
		x = lon.denormalize(Z3(z).d0);
		y = lat.denormalize(Z3(z).d1);
		t = time.denormalize(Z3(z).d2);
	}

	virtual  SeqIndexRange ranges(SeqBoxs seqboxs, SeqTimes seqtimes, int precision, int maxRanges)
	{
		ArrayZRange az;
		for (int i = 0; i < seqboxs.size(); i++)
		{
			ull min = index(seqboxs[i][0], seqboxs[i][1], seqtimes[i][0]);
			ull max = index(seqboxs[i][2], seqboxs[i][3], seqtimes[i][1]);
			az.emplace_back(min, max);
		}

		Z3 zn(3);
		return zn.zranges(az, precision, maxRanges, 7);
	}
};
#pragma endregion


#pragma region XZ3Ordering =XZ3SFC.scala
class XZ3SpaceTimeFillingCurve :public  SpaceTimeFillingCurve 
{
	// 通过 SpaceTimeFillingCurve 继承
	virtual ull index(double x, double y, ull t, bool lenient = false) override
	{
		return ull();
	}
	virtual void invert(ull z, double & x, double & y, ull & t) override
	{
	}
	virtual SeqIndexRange ranges(SeqBoxs seqboxs, SeqTimes seqtimes, int precision, int maxRanges) override
	{
		return SeqIndexRange();
	}
};

#pragma endregion


#pragma region S2Ordering = S2SFC.scala


struct Hilbert {
	//旋转象限
	static void rot(long long  n, long long * x, long long * y, long long  rx, long long  ry) {
		if (ry == 0) {
			if (rx == 1) {
				*x = n - 1 - *x;
				*y = n - 1 - *y;
			}

			//Swap x and y
			long long  t = *x;
			*x = *y;
			*y = t;
		}
	}


	static long long  xy2d(long long  n, long long  x, long long  y) {
		long long  rx, ry, s, d = 0;
		for (s = n / 2; s > 0; s /= 2) {

			//由于s是2的某次方，即在该位为1，而其他位数全为0，
			//则如果x<s，x&s=0，x>=s ,x&s>0
			rx = (x & s) > 0;
			ry = (y & s) > 0;

			//(3 * rx) ^ ry这一块象限排序规则与之前的python代码完全相同
			d += s * s * ((3 * rx) ^ ry);
			rot(s, &x, &y, rx, ry);
		}
		return d;
	}

	//convert d to (x,y)
	static void d2xy(long long  n, long long d, long long * x, long long * y) {
		long long rx, ry, s, t = d;
		*x = *y = 0;
		for (s = 1; s < n; s *= 2) {
			rx = 1 & (t / 2);
			ry = 1 & (t ^ rx);
			rot(s, x, y, rx, ry);
			*x += s * rx;
			*y += s * ry;
			t /= 4;
		}
	}
};


class S2SpaceTimeFillingCurve :public  SpaceTimeFillingCurve
{
	double LonMin = -180;
	double LonMax = 180;
	double LatMin = -90;
	double LatMax = 90;
	S2SpaceTimeFillingCurve(int minLevel,int maxLevel ,int levelMod ,int maxCells)
	{

	}
	// 通过 SpaceTimeFillingCurve 继承
	virtual ull index(double x, double y, ull t, bool lenient = false) override
	{
		return ull();
	}
	virtual void invert(ull z, double & x, double & y, ull & t) override
	{
	}
	virtual SeqIndexRange ranges(SeqBoxs seqboxs, SeqTimes seqtimes, int precision, int maxRanges) override
	{
		return SeqIndexRange();
	}
};

#pragma endregion