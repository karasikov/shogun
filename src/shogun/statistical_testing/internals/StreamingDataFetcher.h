/*
 * Restructuring Shogun's statistical hypothesis testing framework.
 * Copyright (C) 2016  Soumyajit De
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http:/www.gnu.org/licenses/>.
 */

#include <memory>
#include <shogun/lib/common.h>
#include <shogun/statistical_testing/internals/DataFetcher.h>

#ifndef STREMING_DATA_FETCHER_H__
#define STREMING_DATA_FETCHER_H__

namespace shogun
{

class CStreamingFeatures;

namespace internal
{

class DataManager;
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class StreamingDataFetcher : public DataFetcher
{
	friend class DataManager;
public:
	StreamingDataFetcher(CStreamingFeatures* samples);
	virtual ~StreamingDataFetcher();
	void set_num_samples(index_t num_samples);

	virtual void shuffle_features();
	virtual void unshuffle_features();

	virtual void use_fold(index_t i);
	virtual void init_active_subset();

	virtual void start();
	virtual CFeatures* next();
	virtual void reset();
	virtual void end();

	virtual index_t get_num_samples() const;
	virtual const char* get_name() const
	{
		return "StreamingDataFetcher";
	}
private:
	std::shared_ptr<CStreamingFeatures> m_samples;
	bool parser_running;
};
#endif // DOXYGEN_SHOULD_SKIP_THIS
}

}
#endif // STREMING_DATA_FETCHER_H__
