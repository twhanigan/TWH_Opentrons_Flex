def next_tip(
        self,
        num_tips: int = 1,
        starting_tip: Optional[Well] = None,
        *,
        nozzle_map: Optional[NozzleMapInterface] = None,
    ) -> Optional[Well]:
        """
        Find the next valid well for pick-up.

        Determines the next valid start tip from which to retrieve the
        specified number of tips. There must be at least `num_tips` sequential
        wells for which all wells have tips, in the same column.

        :param num_tips: target number of sequential tips in the same column
        :type num_tips: int
        :param starting_tip: The :py:class:`.Well` from which to start search.
                for an available tip.
        :type starting_tip: :py:class:`.Well`
        :return: the :py:class:`.Well` meeting the target criteria, or None
        """
        assert num_tips > 0, f"num_tips must be positive integer, but got {num_tips}"

        well_name = self._core.get_next_tip(
            num_tips=num_tips,
            starting_tip=starting_tip._core if starting_tip else None,
            nozzle_map=nozzle_map,
        )
