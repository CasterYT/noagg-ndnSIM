/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2011-2018  Regents of the University of California.
 *
 * This file is part of ndnSIM. See AUTHORS for complete list of ndnSIM authors and
 * contributors.
 *
 * ndnSIM is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * ndnSIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * ndnSIM, e.g., in COPYING.md file.  If not, see <http://www.gnu.org/licenses/>.
 **/

#include "ndn-consumer-INA.hpp"
#include <fstream>
#include <string>

NS_LOG_COMPONENT_DEFINE("ndn.ConsumerINA");

namespace ns3 {
    namespace ndn {

        NS_OBJECT_ENSURE_REGISTERED(ConsumerINA);


        TypeId
        ConsumerINA::GetTypeId()
        {
            static TypeId tid =
                    TypeId("ns3::ndn::ConsumerINA")
                            .SetGroupName("Ndn")
                            .SetParent<Consumer>()
                            .AddConstructor<ConsumerINA>()

                            .AddAttribute("Beta",
                                          "TCP Multiplicative Decrease factor",
                                          DoubleValue(0.5),
                                          MakeDoubleAccessor(&ConsumerINA::m_beta),
                                          MakeDoubleChecker<double>())

                            .AddAttribute("AddRttSuppress",
                                          "Minimum number of RTTs (1 + this factor) between window decreases",
                                          DoubleValue(0.5), // This default value was chosen after manual testing
                                          MakeDoubleAccessor(&ConsumerINA::m_addRttSuppress),
                                          MakeDoubleChecker<double>())

                            .AddAttribute("ReactToCongestionMarks",
                                          "If true, process received congestion marks",
                                          BooleanValue(true),
                                          MakeBooleanAccessor(&ConsumerINA::m_reactToCongestionMarks),
                                          MakeBooleanChecker())

                            .AddAttribute("UseCwa",
                                          "If true, use Conservative Window Adaptation",
                                          BooleanValue(true),
                                          MakeBooleanAccessor(&ConsumerINA::m_useCwa),
                                          MakeBooleanChecker())
                            .AddAttribute("Window", "Initial size of the window", StringValue("1"),
                                          MakeUintegerAccessor(&ConsumerINA::GetWindow, &ConsumerINA::SetWindow),
                                          MakeUintegerChecker<uint32_t>())

                            .AddAttribute("InitialWindowOnTimeout", "Set window to initial value when timeout occurs",
                                          BooleanValue(true),
                                          MakeBooleanAccessor(&ConsumerINA::m_setInitialWindowOnTimeout),
                                          MakeBooleanChecker())
                            .AddTraceSource("WindowTrace",
                                            "Window that controls how many outstanding interests are allowed",
                                            MakeTraceSourceAccessor(&ConsumerINA::m_window),
                                            "ns3::ndn::ConsumerINA::WindowTraceCallback")
                            .AddTraceSource("InFlight", "Current number of outstanding interests",
                                            MakeTraceSourceAccessor(&ConsumerINA::m_inFlight),
                                            "ns3::ndn::ConsumerINA::WindowTraceCallback");

            return tid;
        }

        ConsumerINA::ConsumerINA()
                : m_inFlight(0)
                , m_ssthresh(std::numeric_limits<double>::max())
                , m_highData(0)
                , m_recPoint(0.0)
        {
            // Open and immediately close the file in write mode to clear it
            std::ofstream file1(windowTimeRecorder, std::ios::out);
            if (!file1.is_open()) {
                std::cerr << "Failed to open the file: " << windowTimeRecorder << std::endl;
            }
            file1.close(); // Optional here since file will be closed automatically
        }

        void
        ConsumerINA::SendInterest(shared_ptr<Name> newName)
        {
            m_inFlight++;
            Consumer::SendInterest(newName);
        }


        void
        ConsumerINA::ScheduleNextPacket()
        {
            if (m_window == static_cast<uint32_t>(0)) {
                Simulator::Remove(m_sendEvent);

                NS_LOG_DEBUG("Error! Window becomes 0!!!!!!");
                //NS_LOG_DEBUG("New event in " << (std::min<double>(0.5, m_rtt->RetransmitTimeout().ToDouble(Time::S))) << " sec");
                m_sendEvent = Simulator::Schedule(Seconds(std::min<double>(0.5, (m_retxTimer*6).GetSeconds())),
                                                  &Consumer::SendPacket, this);
            }
            else if (m_inFlight >= m_window) {
                // do nothing
            }
            else {
                if (m_sendEvent.IsRunning()) {
                    Simulator::Remove(m_sendEvent);
                }
                NS_LOG_DEBUG("Window: " << m_window << ", InFlight: " << m_inFlight);
                m_sendEvent = Simulator::ScheduleNow(&Consumer::SendPacket, this);
            }
        }


        void
        ConsumerINA::StartApplication()
        {
            Consumer::StartApplication();
        }

        void
        ConsumerINA::OnData(shared_ptr<const Data> data)
        {
            Consumer::OnData(data);

            WindowMeasure();

            uint64_t sequenceNum = data->getName().get(-1).toSequenceNumber();

            // Set highest received Data to sequence number
            if (m_highData < sequenceNum) {
                m_highData = sequenceNum;
            }

            if (data->getCongestionMark() > 0) {
                if (m_reactToCongestionMarks) {
                    NS_LOG_DEBUG("Received congestion mark: " << data->getCongestionMark());
                    WindowDecrease();
                }
                else {
                    NS_LOG_DEBUG("Ignored received congestion mark: " << data->getCongestionMark());
                }
            }
            else {
                WindowIncrease();
            }

            if (m_inFlight > static_cast<uint32_t>(0)) {
                m_inFlight--;
            }

            NS_LOG_DEBUG("Window: " << m_window << ", InFlight: " << m_inFlight);

            ScheduleNextPacket();
        }

        void
        ConsumerINA::OnTimeout(std::string nameString)
        {
            WindowDecrease();

            if (m_inFlight > static_cast<uint32_t>(0)) {
                m_inFlight--;
            }

            NS_LOG_DEBUG("Window: " << m_window << ", InFlight: " << m_inFlight);

            Consumer::OnTimeout(nameString);
        }

        void
        ConsumerINA::SetWindow(uint32_t window)
        {
            m_initialWindow = window;
            m_window = m_initialWindow;
        }

        uint32_t
        ConsumerINA::GetWindow() const
        {
            return m_initialWindow;
        }


        void
        ConsumerINA::WindowIncrease()
        {
            if (m_window < m_ssthresh) {
                m_window += 1.0;
            }else {
                m_window += (1.0 / m_window);
            }
            NS_LOG_DEBUG("Window size increased to " << m_window);
        }


        void
        ConsumerINA::WindowDecrease()
        {
            const double diff = m_seq - m_highData;
            BOOST_ASSERT(diff >= 0);

            m_recPoint = m_seq + (m_addRttSuppress * diff);

            // AIMD
            m_ssthresh = m_window * m_beta;
            m_window = m_ssthresh;

            // Window size can't be reduced below initial size
            if (m_window < m_initialWindow) {
                m_window = m_initialWindow;
            }

            NS_LOG_DEBUG("Window size decreased to " << m_window);
        }


        void
        ConsumerINA::WindowMeasure()
        {
            // Open file; on first call, truncate it to delete old content
            std::ofstream file(windowTimeRecorder, std::ios::app);

            if (file.is_open()) {
                file << ns3::Simulator::Now().GetMilliSeconds() << " " << m_window << "\n";  // Write text followed by a newline
                file.close();          // Close the file after writing
            } else {
                std::cerr << "Unable to open file: " << windowTimeRecorder << std::endl;
            }
            //windowMonitor = Simulator::Schedule(MilliSeconds(5), &ConsumerINA::WindowMeasure, this);
        }

    } // namespace ndn
} // namespace ns3
