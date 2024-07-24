//
// Created by 李怿曈 on 16/4/2024.
//

#ifndef NDN_AGGREGATOR_HPP
#define NDN_AGGREGATOR_HPP

#include "ns3/ndnSIM/model/ndn-common.hpp"

#include "ndn-app.hpp"
#include "ModelData.hpp"

#include "ns3/random-variable-stream.h"
#include "ns3/nstime.h"
#include "ns3/data-rate.h"
#include "ns3/traced-value.h"

#include "ns3/ndnSIM/model/ndn-common.hpp"
#include "ns3/ndnSIM/utils/ndn-rtt-estimator.hpp"
#include "ns3/ptr.h"

#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <numeric>
#include <iostream>
#include <sstream>
#include <queue>
#include <utility>
#include <tuple>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/tag.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

namespace ns3{
namespace ndn{

class Aggregator : public App
{
public:

    // New design to get child info's mapping
    void
    AddLinkInfo(const std::string& filename);

    // Calculate response time
    void
    ResponseTimeSum (int64_t response_time);

    int64_t
    GetResponseTimeAverage();

    // Calculate aggregate time
    void
    AggregateTimeSum (int64_t response_time);

    int64_t
    GetAggregateTimeAverage();

    static TypeId GetTypeId();


    Aggregator();
    virtual ~Aggregator(){};


    // Parse the received aggregation tree, important for testing!
    std::pair<std::string, std::set<std::string>> aggTreeProcessSingleString(const std::string& input);

    std::map<std::string, std::set<std::string>> aggTreeProcessStrings(const std::vector<std::string>& inputs);








    // Core part
    virtual void OnInterest(shared_ptr<const Interest> interest);

    virtual void OnData(shared_ptr<const Data> data);

    virtual void ScheduleNextPacket();

    void SendPacket();

    void SendInterest(shared_ptr<Name> newName, uint32_t seq);

    void SetWindow(uint32_t payload);

    uint32_t GetWindow() const;

    void SetSeqMax(uint32_t seqMax);

    uint32_t GetSeqMax() const;

    void WindowIncrease();

    void WindowDecrease();








    void RetransmitInterest(shared_ptr<Name> newName, uint32_t seqNum);

    virtual void OnNack(shared_ptr<const lp::Nack> nack);

    virtual void OnTimeout(uint32_t sequenceNumber);

    virtual void WillSendOutInterest(uint32_t sequenceNumber);

    void aggregate(const ModelData& data, const std::string& dataName);

    ModelData getMean(const std::string& dataName);

    // Override the function in App class to return leaf nodes
    std::map<std::string, std::set<std::string>> getLeafNodes(
            const std::string& key,
            const std::map<std::string, std::vector<std::string>>& treeMap);


protected:
    virtual void StartApplication() override;
    virtual void StopApplication() override;

    void CheckRetxTimeout();

    void SetRetxTimer(Time retxTimer);

    Time GetRetxTimer() const;



protected:
    Ptr<UniformRandomVariable> m_rand;
    Name m_prefix;
    Name m_nexthop;
    Name m_nexttype;
    Time m_interestLifeTime;
    std::string m_nodeprefix;


public:
    typedef void (*LastRetransmittedInterestDataDelayCallback)(Ptr<App> app, uint32_t seqno, Time delay, int32_t hopCount);
    typedef void (*FirstInterestDataDelayCallback)(Ptr<App> app, uint32_t seqno, Time delay, uint32_t retxCount, int32_t hopCount);

protected:
    // Important queue definition, tuple: interest name(string), seq(uint32_t), resKey(string)
    std::queue<std::tuple<shared_ptr<Name>, uint32_t, std::string>> interestQueue;

    // Window design, important!!!!!!
    uint32_t m_initialWindow;
    TracedValue<double> m_window;;
    TracedValue<uint32_t> m_inFlight;
    bool m_setInitialWindowOnTimeout;

    // AIMD design, important!!!!!!
    double m_ssthresh;
    bool m_useCwa;
    uint32_t m_highData;
    double m_recPoint;
    double m_beta;
    double m_addRttSuppress;
    bool m_reactToCongestionMarks;


    uint32_t m_seq;      ///< @brief currently requested sequence number
    uint32_t m_seqMax;   ///< @brief maximum number of sequence number
    EventId m_sendEvent; ///< @brief EventId of pending "send packet" event
    Time m_retxTimer;    ///< @brief Currently estimated retransmission timer
    EventId m_retxEvent; ///< @brief Event to check whether or not retransmission should be performed

    Ptr<RttEstimator> m_rtt; ///< @brief RTT estimator

    Time m_offTime;          ///< \brief Time interval between packets
    Name m_interestName;     ///< \brief NDN Name of the Interest (use Name)
    Time m_freshness;
    uint32_t m_signature;
    Name m_keyLocator;

/*    // Designed for timeout and retransmission, need to re-design if name definition changes
    std::map<uint32_t, std::string> m_newSeq_Name;

    // Designed for data aggregation, currently only for binary tree
    std::map<uint32_t, std::vector<std::string>> map_agg_newSeq_Name;
    std::map<std::string, std::vector<std::string>> map_agg_oldSeq_newName;
    std::map<std::string, std::vector<std::string>> map_agg_final;*/

    // Retransmission
    std::map<uint32_t, std::vector<std::string>> m_newSeq_Name;

    // Check whether aggregation has finished
    //std::map<std::string, std::vector<std::string>> map_agg_final;
    std::map<std::string, std::vector<std::string>> map_agg_oldSeq_newName;
    std::map<uint32_t, std::string> m_agg_newDataName;

    // Designed for actual aggregation operations
    std::map<std::string, std::vector<float>> sumParameters; // result after performing mean average aggregation
    std::map<std::string, int> count; // count of aggregation

    // Check Timeout, important!!!!!!!!!!!!!
    std::map<uint32_t, ns3::Time> m_timeoutCheck;
    ns3::Time m_timeoutThreshold;


    // Designed for response time measurement
/*    std::map<uint32_t, ns3::Time> currentTime;
    std::map<std::string, ns3::Time> aggregateStartTime;
    std::map<uint32_t, ns3::Time> responseTime;
    std::map<std::string, ns3::Time> aggregateTime;*/

    std::map<std::string, ns3::Time> currentTime;
    std::map<uint32_t, ns3::Time> aggregateStartTime;
    std::map<std::string, ns3::Time> responseTime;
    std::map<uint32_t, ns3::Time> aggregateTime;
    int64_t totalAggregateTime;
    int64_t totalResponseTime;
    int round;
    int iteration;

    // Topology file name, important!!!!!!!!!!
    std::string fileAgg = "scratch/agg-tree-tier5.txt";
    std::string fileAggBottleneck = "scratch/agg-tree-tier5-bottleneck.txt";


    // Store aggregation logic, the most important!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std::map<std::string, std::set<std::string>> aggregationMap;


    struct RetxSeqsContainer : public std::set<uint32_t> {
    };

    RetxSeqsContainer m_retxSeqs; ///< \brief ordered set of sequence numbers to be retransmitted


    struct SeqTimeout {
        SeqTimeout(uint32_t _seq, Time _time)
            : seq(_seq)
            , time(_time)
        {
        }

        uint32_t seq;
        Time time;
        };

class i_seq {};

class i_timestamp {};


struct SeqTimeoutsContainer
                : public boost::multi_index::
                multi_index_container<SeqTimeout,
                        boost::multi_index::
                        indexed_by<boost::multi_index::
                        ordered_unique<boost::multi_index::tag<i_seq>,
                  boost::multi_index::
                  member<SeqTimeout, uint32_t,
                          &SeqTimeout::seq>>,
        boost::multi_index::
        ordered_non_unique<boost::multi_index::
        tag<i_timestamp>,
        boost::multi_index::
        member<SeqTimeout, Time,
                &SeqTimeout::time>>>> {
        };

SeqTimeoutsContainer m_seqTimeouts; ///< \brief multi-index for the set of SeqTimeout structs
SeqTimeoutsContainer m_seqLastDelay;
SeqTimeoutsContainer m_seqFullDelay;
std::map<uint32_t, uint32_t> m_seqRetxCounts;


TracedCallback<Ptr<App> /* app */, uint32_t /* seqno */, Time /* delay */, int32_t /*hop count*/> m_lastRetransmittedInterestDataDelay;
TracedCallback<Ptr<App> /* app */, uint32_t /* seqno */, Time /* delay */,uint32_t /*retx count*/, int32_t /*hop count*/> m_firstInterestDataDelay;

};

} // namespace ndn
} // namespace ns3


#endif //APPS_NDN_AGGREGATOR_HPP
