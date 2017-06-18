import heapq


class Item(object):
    def __init__(self, name, weight):
        """Item.
        Args:
            name (str): Name of this item.
            weight (int): Weight of this item.
            values (list): list of value parameters of this item.
        """

        self._name = name
        self._weight = weight

    @property
    def name(self):
        """Return name of Item.
        Return:
            str: Name of the Item.
        """

        return self._name

    @property
    def weight(self):
        """Return weight of the Item.
        Return:
            int: Weight of the Item.
        """

        return self._weight

    def __lt__(self, other):
        """Overload magic method.
        Return:
            bool: This weight is smaller than other weight.
        """

        return self._weight < other._weight

    def __cmp__(self, other):
        """Overload magic method cmp
        Return:
            bool: This weight is smaller than other weight.
        """

        return self._weight < other._weight

    def __eq__(self, other):
        """Overload magic method eq
        Return:
            bool: This weight is equal to other weight
        """

        return self._weight == other._weight


class Bin(object):
    def __init__(self, capacity):
        """Initialize Bin instance
        """

        self._items = []
        self._utilization = 0
        self._capacity = capacity

    @property
    def name(self):
        """Return name.
        Return:
            str: Return name of the Bin.
        """

        return self._name

    @name.setter
    def name(self, name):
        """Set name.
        Args:
            str: name of the Bin.
        """

        self._name = name

    @property
    def capacity(self):
        """Return capacity
        """

        return self._capacity

    @property
    def total_weight(self):
        """Return subtotal of the weights
        """

        return sum([i.weight for i in self._items])

    @property
    def utilization(self):
        """Return utilization.
        Return:
            float: Return utilization in percentage.
        """

        total_weight = sum([i.weight for i in self._items])
        return round((total_weight / self._capacity) * 100, 2)


    def push(self, item):
        """Push Item into heap.
        Args:
            item (Item): generic item.
        """

        heapq.heappush(self._items, item)

    def pop(self):
        """Pop Item from the heap
        Return:
            Item: smallest generic item.
        """

        return heapq.heappop(self._items)

    def remove(self, name):
        """Remove an Item by name from the heap.
        Args:
            name (str): name of the item.
        """

        index = 0
        for i, item in enumerate(self._items):
            index = i
            if item.name == name:
                break
        del self._items[index]
        heapq.heapify(self._items)


    def get_items(self, attribute=None, value=None):
        """Search for item in the bin by attribute
        If attribute or value is not specified, this will return all
        Items in the Bin.
        Args:
            attribute (str)[optional]: one of the attributes on the Ttem.
            value (mixed)[optional]: value of the attribute.
        Returns:
            list: list of Item(s).
        """
        if attribute and value:
            return [i for i in self._items if getattr(i, attribute) == value]
        return self._items


class Binpacker(object):
    def __init__(self, capacity):
        """Create Binpacker instance.
        Args:
            capacity (int): Capacity of one container.
        """

        self._capacity = capacity
        self._bins = []
        self._items = []

    @property
    def items(self):
        """Return items.
        Return:
            list: Return list of Item objects.
        """

        return self._items

    @items.setter
    def items(self, items):
        """Set items.
        Args:
            list: list of Item objects.
        """

        for i in items:
            if i.weight > self._capacity:
                raise Exception(
                    'Weight of Item: "{}" exceeds the maximum '
                    'bin capacity: {}.'.format(
                        i.name, self._capacity))
        self._items = items

    @property
    def bins(self):
        """Return Bin objects.
        """

        return self._bins

    @bins.setter
    def bins(self, bins):
        """Set Bins.
        Args:
            list: list of Bin.
        """

        self._bins = bins

    def get_truth_table(self, capacity, items):
        """Generate truth table.
        """

        m = [
            [True for j in range(capacity + 1)]
            for i in range(len(items))
        ]
        for index, item in enumerate(items, 1):
            i = index - 1
            for j in range(0, capacity + 1):
                if not i:
                    if j != item.weight and j > 0:
                        m[i][j] = False
                else:
                    if j < item.weight:
                        m[i][j] = m[i - 1][j]
                    else:
                        m[i][j] = m[i - 1][j] or m[i - 1][j - item.weight]
        return m

    def _pick_items(self, truth_table):
        """Pick the right Items to be packed into a Bin
        Args:
            truth_table (list): current truth table from best fit algorithm.
        Returns:
            list: list of indices of the picked Items.
        """

        k = len(truth_table) - 1
        picked_items_indices = []
        if k >= 0:
            # Get the heaviest subtotal weight (best utilization)
            # a bin can be loaded.
            j = max([index for index, x in enumerate(truth_table[k]) if x])

            while k >= 0:
                if not k:
                    if j > 0:
                        picked_items_indices.append(k)
                else:
                    if not truth_table[k - 1][j]:
                        picked_items_indices.append(k)
                        j -= self._items[k].weight
                k -= 1
        return picked_items_indices

    def _move_items_to_bin(self, list_of_items_indices, bin_index):
        """Move Items into Bin
        Args:
            list_of_items_indices (list): list of indices of the Items.
            bin_index (int): index of the destination Bin.
        """

        for i in list_of_items_indices:
            self._bins[bin_index].push(self._items[i])
            del self._items[i]

    def pack_items(self,k):
        """Pack items into bin(s).
        Calculate minimum number of bin(s) needed and what items
        go to which bin. At the end of this method, all items
        will be packed into bin(s). The state is going to be updated
        on demand, meaning each Bin can be in the middle of being emptied
        out when it is being packed. Therefore,
        """

        self._items = sorted(self._items)
        # First, go through existing bin to see if they can
        # accommodate new weights.
        for index, old_bin in enumerate(self._bins):
            if old_bin.utilization == 100.00:
                continue
            remaining_space = old_bin.capacity - old_bin.total_weight
            m = self.get_truth_table(remaining_space, self._items)
            picked_items = self._pick_items(m)
            self._move_items_to_bin(picked_items, index)

        # If there is still remaining Items, we need to get a new Bin.
        while len(self._items) > 0:
            m = self.get_truth_table(self._capacity, self._items)
            new_bin = Bin(self._capacity)
            new_bin.name = '[NEW BIN {}]'.format(len(self._bins))
            self._bins.append(new_bin)
            if (len(self._bins) > k):
                return False
            else:
                bin_index = len(self._bins) - 1
                picked_items = self._pick_items(m)
                self._move_items_to_bin(picked_items, bin_index)
        if(len(self._items)== 0):
            return True